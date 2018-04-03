package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::PrePipelineChecks;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

my $previous_version_species_suffix = '_previous_version';

sub run {
    my $self = shift;

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    my @error_msg;
    
    my $species = $self->param('species');
    
    use Bio::EnsEMBL::Registry;
    my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    
    eval {
      my $previous_version_funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( 
        $species . $previous_version_species_suffix, 'funcgen' 
      );
    };
    if ($@) {
      push @error_msg, "Can't get previous version database for $species!";
    }
    
    my $analysis_adaptor = $funcgen_db_adaptor->get_AnalysisAdaptor;
    
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
      push @error_msg, "The following analyses are required by the pipeline, but missing in the database: " . join ', ', @missing_analyses;
    }
    
    my @logic_names_used_in_filenames = @required_logic_name;
    my @analyses_with_problem_names;
    
    foreach my $logic_name (@logic_names_used_in_filenames) {
    
      my $has_spaces = $logic_name =~ / /;
      
      if ($has_spaces) {
        push @analyses_with_problem_names, "The analysis with logic_name '$logic_name' has spaces in it.";
      }
    }
    if (@analyses_with_problem_names) {
      push @error_msg, "The following analyses have logic names that will lead to issues, because they are used in file names on the ftp site:\n" 
        . join "\n", map { "    - " . $_ } @analyses_with_problem_names;
    }

    # Find read files that have been registered, but don't exist on the file 
    # system.
    #
    my $read_file_adaptor = $funcgen_db_adaptor->get_ReadFileAdaptor;
    
    my $registered_file = [ map { $_->file } @{$read_file_adaptor->fetch_all} ];
    
    foreach my $current_file (@$registered_file) {
      if (! -e $current_file) {
        push @error_msg, "Can't find read file: $current_file" 
      }
    }
    
    my $dbc = $funcgen_db_adaptor->dbc;
    
    # Find duplicate feature types
    #
    my $find_duplicate_feature_types = 'select name, max(feature_type_id), group_concat(feature_type_id), count(feature_type_id) c from feature_type where name != "WCE" group by name having c>1;';
    
    my $sth = $dbc->prepare($find_duplicate_feature_types);
    $sth->execute;
    my $duplicate_feature_types = $sth->fetchall_hashref('name');
    my @duplicate_feature_type_names = keys %$duplicate_feature_types;
    
    if (@duplicate_feature_type_names) {
      push @error_msg, "The following feature types are referenced by experiments, but they are not unique: " . join ', ', @duplicate_feature_type_names;
    }

    # Find feature types without SO accessions
    #
    my $find_no_SO_feature_types = 'select feature_type.name from experiment join feature_type using (feature_type_id) where so_accession is null and feature_type.name!="WCE"';
    
    my $sth = $dbc->prepare($find_no_SO_feature_types);
    $sth->execute;
    my $feature_types_no_SO = $sth->fetchall_arrayref;
    
    if (@$feature_types_no_SO) {
      push @error_msg, "The following feature types have no sequence ontology accession: " . join ', ', @$feature_types_no_SO;
    }
    
    # Find duplicate configurations within the same experiment
    #
    my $find_duplicate_configurations = 'select group_concat(read_file_experimental_configuration_id) as invalid_configuration, count(read_file_experimental_configuration_id) c from read_file_experimental_configuration group by experiment_id, biological_replicate, technical_replicate, multiple, paired_end_tag having c>1;';
    
    my $sth = $dbc->prepare($find_duplicate_configurations);
    $sth->execute;
    my $invalid_configurations = $sth->fetchall_hashref('invalid_configuration');
    
    if (keys %$invalid_configurations) {
      push @error_msg, "The following read file configuration ids are duplicates within the same experiment: " . join ', ', map { '('.$_.')' } keys %$invalid_configurations;
    }
    
    # Feature types have SO accessions that are in the ontology database
    #
    my $find_feature_type_ids_needing_valid_SO = 'select distinct so_accession from experiment join feature_type using (feature_type_id) where so_accession is not null and feature_type.name!="WCE"';
    
    my $sth = $dbc->prepare($find_feature_type_ids_needing_valid_SO);
    $sth->execute;
    my $so_accessions = $sth->fetchall_arrayref;
    
    my $ontology_term_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
    
    my @not_found;
    foreach my $so_accession (@$so_accessions) {
      $so_accession = $so_accession->[0];
      my $ontology_term = $ontology_term_adaptor->fetch_by_accession($so_accession);
      
      if (! defined $ontology_term) {
        push @not_found, $so_accession
      }
    }
    if (@not_found) {
      push @error_msg, "The following sequence ontology accessions can't be fetched from the ontology database: " . join ', ', @not_found;
    }
    
    my $find_signal_control_epigenome_mismatches = 'select * from experiment s join experiment c on (s.control_id=c.experiment_id) and s.epigenome_id != c.epigenome_id';
    
    my $sth = $dbc->prepare($find_signal_control_epigenome_mismatches);
    $sth->execute;
    my $signal_control_epigenome_mismatches = $sth->fetchall_arrayref;
    if (@$signal_control_epigenome_mismatches) {
        push @error_msg, "Signal control mismatches!";
    }
    
    # Check for orphans
    #
    my $count_orphans_sql = 'select count(*) as c from alignment_read_file left join read_file using (read_file_id) where read_file.read_file_id is null';
    $sth = $dbc->prepare($count_orphans_sql);
    $sth->execute;
    
    my $x = $sth->fetchall_arrayref;
    my $num_orphans = $x->[0]->[0];
    
    if ($num_orphans) {
        push @error_msg, "There are $num_orphans orphan entries in the alignment_read_file table."; 
    }

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
    );
    
    foreach my $current_program (@programs_expected_in_path) {
      system("which $current_program > /dev/null");
      if ($?) {
    push @error_msg, "Can't find $current_program in path.";
      }
    }
    
    system("which Rscript > /dev/null");    
    if ($?) {
      push @error_msg, "Can't find Rscript in path.";
    } else {
    
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
      
      if ($version_string_found) {
      
    my $version = $1;
    
    # If version string starts with 3.2, we are probably ok.
    #
    my $version_ok =
         $version =~ /^3\.2/
      || $version =~ /^3\.3/
      || $version =~ /^3\.4/
    ;
    
    if (!$version_ok) {
      push @error_msg, "Rscript on PATH has version $version, this may cause trouble.";
    }   
      } else {
    push @error_msg, "Can't find version from Rscript in string: '$version_string'";
      }
    }
    if (@error_msg) {
    
      my $err_string = join "\n", map { '  - ' . $_ } @error_msg;

      my $complete_error_message 
        = 
        "\n-------------------------------------------------------------------------------\n"
        . "Prepipeline checks have failed with the following errors:\n\n"
        . $err_string
        . "\n\nPlease fix these first before continuing with the pipeline."
        . "\n-------------------------------------------------------------------------------\n";
        
      $self->throw($complete_error_message);
    }
    $logger->info("Pre pipeline checks have completed successfully. You can now run the rest of the pipeline.\n");
}

1;
