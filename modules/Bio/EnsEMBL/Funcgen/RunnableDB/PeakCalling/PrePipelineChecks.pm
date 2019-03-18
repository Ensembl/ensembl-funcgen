package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::PrePipelineChecks;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

my $previous_version_species_suffix = '_previous_version';

sub species {
  my $self = shift;
  return $self->param('species');
}

sub funcgen_adaptor {
  my $self = shift;
  return Bio::EnsEMBL::Registry->get_DBAdaptor(
      $self->species,
      'funcgen'
   );
}

sub read_file_adaptor {
  my $self = shift;
  return $self->funcgen_adaptor->get_ReadFileAdaptor;
}

sub funcgen_dbc {
  my $self = shift;
  return $self->funcgen_adaptor->dbc;
}

sub run {
    my $self = shift;

    my @error_msg;
    
    my $dbc = $self->funcgen_dbc;
    
    my @error_messages = (
      $self->check_analyses,
      $self->check_registry_has_previous_version_db,
      $self->check_read_files_exist,
      $self->check_read_file_names_ok,
      $self->check_unused_controls($dbc),
      #$self->check_feature_types_unique($dbc),
      #$self->check_feature_types_have_so_accessions($dbc),
      $self->check_read_file_configurations_unique_per_experiment($dbc),
      $self->check_so_accessions_in_ontology_db($dbc),
      #$self->check_signal_control_mismatches($dbc),
      #$self->check_alignment_read_file_orphans($dbc),
      #$self->check_duplicate_read_file_names($dbc),
      $self->check_environment_variables_exported,
      $self->check_R,
      $self->check_programs_in_path,
      #$self->check_unused_control_reads_files_although_they_do_control_something($dbc),
    );

    push @error_msg, grep { defined $_ } @error_messages;
    
    if (@error_msg) {
    
      my $err_string = join "\n", map { "\n-------------------\n" . $_ } @error_msg;

      my $complete_error_message 
        = 
        "\n-------------------------------------------------------------------------------\n"
        . "Prepipeline checks have failed with the following errors:\n\n"
        . $err_string
        . "\n\nPlease fix these first before continuing with the pipeline."
        . "\n-------------------------------------------------------------------------------\n";
        
      $self->throw($complete_error_message);
    }
    $self->say_with_header("Pre pipeline checks have completed successfully. You can now run the rest of the pipeline.", 1);
}

sub check_R {
  my $self = shift;

  my @error_msg;

  system("which Rscript > /dev/null");
  if ($?) {
    return "Can't find Rscript in path.";
  }
    
  # Should be something like 
  # "R scripting front-end version 3.2.2 (2015-08-14)"
  #
  # The one in here
  #
  # /software/R-3.2.2/bin/
  #
  # should work. To make life more interesting, Rscript outputs the 
  # version string to stderr.
  #
  my $version_string = `Rscript --version 2>&1`;
  
  my $version_string_found = $version_string =~ /R scripting front-end version ([^ ]+?) /;
  
  if (! $version_string_found) {
    return "Can't find version from Rscript in string: '$version_string'";
  }

  my $version = $1;
  # If version string starts with 3.2, we are probably ok.
  #
  my $version_ok =
       $version =~ /^3\.2/
    || $version =~ /^3\.3/
    || $version =~ /^3\.4/
  ;

  if (! $version_ok) {
    return "Rscript on PATH has version $version, this may cause trouble.";
  }
  return;
}

sub check_programs_in_path {
  my $self = shift;

  my @error_msg;

  my @programs_expected_in_path = qw(
    bwa
    R
    samtools
    run_spp.R
    load_phantom_peak_file.pl
    load_argenrich_qc_file.pl
    load_fastqc_summary_file.pl
    load_samtools_flagstats.pl
    argenrichformregions.pl
    argenrich_with_labels_and_rerunnable.R
    idr
    populate_meta_coord.pl
    export_regulatory_features_to_bed.pl
  );

  foreach my $current_program (@programs_expected_in_path) {
    system("which $current_program > /dev/null");
    if ($?) {
      push @error_msg, "Can't find the script \"$current_program\" in the search path, but it is used in the regulatory build pipeline.";
    }
  }
  
  my $error_msg;
  if (@error_msg) {
    $error_msg = join "\n", @error_msg;
  }
  
  return $error_msg;
}

sub check_environment_variables_exported {
  my $self = shift;

  my @error_msg;

  my @exported_variables = qw( PERL5LIB R_LIBS );
  my $cmd;
  
  ENVIRONMENT_VARIABLE:
  foreach my $exported_variable (@exported_variables) {
  
    $cmd = qq(set | grep -v BASH_EXECUTION_STRING | grep $exported_variable);
    my $set_command = `$cmd`;
    chomp $set_command;
    
    if (! $set_command) {
      push @error_msg, "$exported_variable has not been set!";
      next ENVIRONMENT_VARIABLE;
    }

    #
    # If "foo" has been exported, the command
    #     
    # export | grep foo
    #
    # will yield this:
    #
    # declare -x foo="bar"
    # 
    $cmd = qq(export | grep $exported_variable);
    my $declare_command = `$cmd`;
    
    if (! $declare_command) {
      push @error_msg, "$exported_variable has not been exported!";
    }
  }

  my $error_msg;
  if (@error_msg) {
    $error_msg = join "\n", $error_msg;
  }
  return $error_msg;
}

sub check_analyses {
  my $self = shift;

  my @error_msg;

  my $analysis_adaptor = $self->funcgen_adaptor->get_AnalysisAdaptor;
  
  my @required_logic_name = (
    ENSEMBL_SINGLE_END_ALIGNMENT_ANALYSIS,
    ENSEMBL_PAIRED_END_ALIGNMENT_ANALYSIS,
    ENSEMBL_HODGEPODGE_ALIGNMENT_ANALYSIS,
    ENSEMBL_REMOVE_DUPLICATES_ANALYSIS,
    ENSEMBL_BROAD_PEAK_CALLING_ANALYSIS,
    ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_PERMISSIVE,
    ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_DEFAULT
  );
  my @missing_analyses;
  foreach my $logic_name (@required_logic_name) {
  
    my $required_analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);
    
    if (! defined $required_analysis) {
      push @missing_analyses, $logic_name;
    }
  }
  if (@missing_analyses) {
    push 
      @error_msg, 
      "The following analyses are required by the pipeline, but missing "
      . "in the database: " 
      . join ', ', @missing_analyses;
  }
  
  my @logic_names_used_in_filenames = @required_logic_name;
  my @analyses_with_problem_names;
  
  foreach my $logic_name (@logic_names_used_in_filenames) {
  
    my $has_spaces = $logic_name =~ / /;
    
    if ($has_spaces) {
      push 
        @analyses_with_problem_names, 
        "The analysis with logic_name '$logic_name' has spaces in it.";
    }
  }
  if (@analyses_with_problem_names) {
    push @error_msg, "The following analyses have logic names that will lead "
      . "to issues, because they are used in file names on the ftp site:\n" 
      . $self->bulletpointify(\@analyses_with_problem_names)
    ;
  }
  
  my $error_msg;
  if (@error_msg) {
    $error_msg = join "\n", $error_msg;
  }
  return $error_msg;
}

sub check_registry_has_previous_version_db {
  my $self = shift;

  my $error_msg;
  my $species = $self->species;
  
  eval {
    my $previous_version_funcgen_adaptor 
      = Bio::EnsEMBL::Registry->get_DBAdaptor( 
          $species . $previous_version_species_suffix, 
          'funcgen' 
        );
  };
  if ($@) {
    $error_msg = "Can't get previous version database for $species!";
  }
  return $error_msg;
}

sub check_unused_controls {
  my $self = shift;
  my $dbc  = shift;
  
  my $error_msg;
  
  my $sql = qq{
    select 
      control.name
    from 
      experiment control 
      left join experiment sig on (
          sig.control_id = control.experiment_id 
      ) 
    where 
      control.is_control = 1 
      and sig.experiment_id is null
   order by 
      control.name
  };
  my $sth = $dbc->prepare($sql);
  $sth->execute;
  my $x = $sth->fetchall_hashref('name');
  
  my @unused_control_names = sort alphabetically keys %$x;
  
  if (@unused_control_names) {
    $error_msg 
      = "The following " . scalar @unused_control_names . " control experiments are never used:\n" 
        . $self->bulletpointify(\@unused_control_names);
  }
  return $error_msg;
}

sub check_read_files_exist {

  my $self = shift;

  # Find read files that have been registered, but don't exist on the file 
  # system.
  #
  my $read_file_adaptor = $self->read_file_adaptor;
  
  my @all_read_files = @{$read_file_adaptor->fetch_all};
  
  my $registered_file = [ map { $_->file } @all_read_files ];
  
  my @error_msg;
  
  foreach my $current_file (@$registered_file) {
    if (! -e $current_file) {
        push @error_msg, "Can't find read file: $current_file" 
    }
  }
  if (@error_msg) {
    return join "\n", @error_msg;
  }
  return;
}

sub check_read_file_names_ok {
  my $self = shift;
  
  my $read_file_adaptor = $self->read_file_adaptor;
  my @all_read_files = @{$read_file_adaptor->fetch_all};
  
  my @error_msg;
  
  foreach my $read_file (@all_read_files) {
      
      my $is_ok = $read_file->name =~  /^[a-zA-Z0-9_\+\-\:\.]+$/;
      if (! $is_ok) {
          push 
            @error_msg, 
            "Read file name can't be used for creating directories or names: " 
            . $read_file->name 
      }
  }
  if (@error_msg) {
    return join "\n", @error_msg;
  }
  return;
}

# sub check_feature_types_unique {
#   my $self = shift;
#   my $dbc  = shift;
#   
#   # Find duplicate feature types
#   #
#   my $find_duplicate_feature_types = '
#     select 
#         feature_type.name, 
#         max(feature_type_id), 
#         group_concat(feature_type_id), 
#         count(distinct feature_type_id) c 
#     from 
#         feature_type join experiment using (feature_type_id) 
#     where 
#         feature_type.name != "WCE" 
#     group by 
#         feature_type.name 
#     having c>1;
#   ';
#   
#   my $sth = $dbc->prepare($find_duplicate_feature_types);
#   $sth->execute;
#   my $duplicate_feature_types = $sth->fetchall_hashref('name');
#   my @duplicate_feature_type_names = keys %$duplicate_feature_types;
#   
#   my @error_msg;
#   if (@duplicate_feature_type_names) {
#     push 
#       @error_msg, 
#       "The following feature types are referenced by experiments, but they "
#       . "are not unique: " 
#       . join ', ', @duplicate_feature_type_names;
#   }
# 
#   if (@error_msg) {
#     return join "\n", @error_msg;
#   }
#   return;
# }

# sub check_feature_types_have_so_accessions {
#   my $self = shift;
#   my $dbc  = shift;
#   
#   # Find feature types without SO accessions
#   #
#   my $find_no_SO_feature_types = '
#     select 
#       feature_type.name 
#     from 
#       experiment 
#       join feature_type using (feature_type_id) 
#     where 
#       so_accession is null 
#       and feature_type.name != "WCE"
#     '
#   ;
#   
#   my $sth = $dbc->prepare($find_no_SO_feature_types);
#   $sth->execute;
#   my $feature_types_no_SO = $sth->fetchall_arrayref;
#   
#   my @error_msg;
#   if (@$feature_types_no_SO) {
#     push 
#       @error_msg, 
#       "The following feature types have no sequence ontology accession: " 
#       . join ', ', @$feature_types_no_SO;
#   }
# 
#   if (@error_msg) {
#     return join "\n", @error_msg;
#   }
#   return;
# }

sub check_so_accessions_in_ontology_db {
  my $self = shift;
  my $dbc  = shift;
  
  # Feature types have SO accessions that are in the ontology database
  #
  my $find_feature_type_ids_needing_valid_SO = '
    select 
      distinct so_accession 
    from 
      experiment 
      join feature_type using (feature_type_id) 
    where 
      so_accession is not null 
      and feature_type.name!="WCE"
    '
  ;
  
  my $sth = $dbc->prepare($find_feature_type_ids_needing_valid_SO);
  $sth->execute;
  my $so_accessions = $sth->fetchall_arrayref;
  
  my $ontology_term_adaptor 
    = Bio::EnsEMBL::Registry->get_adaptor( 
        'Multi', 
        'Ontology', 
        'OntologyTerm' 
      );
  
  my @not_found;
  foreach my $so_accession (@$so_accessions) {
    $so_accession = $so_accession->[0];
    my $ontology_term 
      = $ontology_term_adaptor->fetch_by_accession($so_accession);
    
    if (! defined $ontology_term) {
      push @not_found, $so_accession
    }
  }
  if (@not_found) {
    return
      "The following sequence ontology accessions can't be "
      . "fetched from the ontology database: " 
      . join ', ', @not_found
    ;
  }
  return;
}

# sub check_signal_control_mismatches {
#   my $self = shift;
#   my $dbc  = shift;
#   
#   my @error_msg;
#   
#   my $find_signal_control_epigenome_mismatches = '
#     select 
#       * 
#     from 
#       experiment s 
#       join experiment c on (
#         s.control_id = c.experiment_id 
#         and s.epigenome_id != c.epigenome_id
#       )
#   ';
#   
#   my $sth = $dbc->prepare($find_signal_control_epigenome_mismatches);
#   $sth->execute;
#   my $signal_control_epigenome_mismatches = $sth->fetchall_arrayref;
#   if (@$signal_control_epigenome_mismatches) {
#     return 
#       "Signal control mismatches! Try:\n\n"
#       . $find_signal_control_epigenome_mismatches;
#   }
#   return;
# }

sub check_read_file_configurations_unique_per_experiment {
  my $self = shift;
  my $dbc  = shift;
  
  # Find duplicate configurations within the same experiment
  #
  my $find_duplicate_configurations = '
    select 
      group_concat(read_file_experimental_configuration_id) as invalid_configuration, 
      count(read_file_experimental_configuration_id) c 
    from 
      read_file_experimental_configuration 
    group by 
      experiment_id, 
      biological_replicate, 
      technical_replicate, 
      multiple, 
      paired_end_tag 
    having 
      c>1
    ;
  ';
  
  my $sth = $dbc->prepare($find_duplicate_configurations);
  $sth->execute;
  my $invalid_configurations = $sth->fetchall_hashref('invalid_configuration');
  
  if (keys %$invalid_configurations) {

    return 
      "The following read file configuration ids are "
      . "duplicates within the same experiment: " 
      . join 
          ', ', 
          map { '('.$_.')' } keys %$invalid_configurations
    ;
  }
  return;
}

# sub check_duplicate_read_file_names {
#   my $self = shift;
#   my $dbc  = shift;
#   
#   my $sql = '
#     select 
#       read_file.name as read_file_name,
#       group_concat(experiment.name) experiment_names,
#       count(read_file_id) count
#     from 
#       experiment 
#           join read_file_experimental_configuration using (experiment_id)
#           join read_file using (read_file_id)
#     group by
#       read_file_name
#     having 
#       count>1
#     ;
#   ';
#   my $sth = $dbc->prepare($sql);
#   $sth->execute;
#   
#   my $hash = $sth->fetchall_hashref('read_file_name');
#   my @duplicate_read_file_names = sort alphabetically keys %$hash;
#   
#  
#   if (@duplicate_read_file_names) {
#       return 
#         "There are ".scalar @duplicate_read_file_names." duplicate read file names in the alignment_read_file table:\n"
#         . $self->bulletpointify(\@duplicate_read_file_names)
#       ;
#   }
#   return;
# }

# sub check_alignment_read_file_orphans {
#   my $self = shift;
#   my $dbc  = shift;
#   
#   my $count_orphans_sql = '
#     select 
#       read_file.name
#     from 
#       read_file left join alignment_read_file using (read_file_id) 
#     where 
#       alignment_read_file.alignment_read_file_id is null
#     order by
#       read_file.name
#   ';
#   my $sth = $dbc->prepare($count_orphans_sql);
#   $sth->execute;
#   
#   my $orphan_read_file_names_hash = $sth->fetchall_hashref('name');
#   my @orphan_read_file_names = sort alphabetically keys %$orphan_read_file_names_hash;
#   
#   if (@orphan_read_file_names) {
#       return 
#         "There are ".scalar @orphan_read_file_names." orphan entries in the alignment_read_file table:\n"
#         . $self->bulletpointify(\@orphan_read_file_names)
#       ;
#   }
#   return;
# }

# sub check_unused_control_reads_files_although_they_do_control_something {
#   my $self = shift;
#   my $dbc  = shift;
#   
#   my $count_orphans_sql = '
#     select 
#       distinct control.name
#     from 
#       experiment control 
#       join experiment sig on (
#           sig.control_id = control.experiment_id 
#       )
#       join read_file_experimental_configuration on (read_file_experimental_configuration.experiment_id = control.experiment_id)
#       join read_file using (read_file_id)
#       left join alignment_read_file using (read_file_id)
#     where 
#       control.is_control = 1
#       and alignment_read_file.alignment_read_file_id is null
#     order by 
#       control.name
#   ';
#   my $sth = $dbc->prepare($count_orphans_sql);
#   $sth->execute;
#   
#   my $orphan_read_file_names_hash = $sth->fetchall_hashref('name');
#   my @orphan_read_file_names = sort alphabetically keys %$orphan_read_file_names_hash;
#   
#   if (@orphan_read_file_names) {
#       return 
#         "There are " . scalar @orphan_read_file_names . " unused control reads files although they do control something:\n"
#         . $self->bulletpointify(\@orphan_read_file_names)
#       ;
#   }
#   return;
# }

sub alphabetically {
  return uc($a) cmp uc($b)
}

sub bulletpointify {
  my $self = shift;
  my $list = shift;
  
  my $bulletpoint = '  - ';
  
  return join "\n", map { $bulletpoint . $_ } @$list;
}

1;
