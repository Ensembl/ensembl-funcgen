package Bio::EnsEMBL::Funcgen::Hive::ErsaPrePipelineChecks;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;

sub run {
    my $self = shift;

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    my @error_msg;
    
    my $out_db_connetion_details = $self->param('out_db');

    use Bio::EnsEMBL::DBSQL::DBAdaptor;
    my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$out_db_connetion_details);
    my $dbc = $dba->dbc;
    
    # Check that the values in the replicate column aren't duplicated.
    
    my $find_non_unique_replicate_specifications = 'select epigenome_id, experiment_id, biological_replicate, technical_replicate, count(*) c from input_subset where is_control=0 group by epigenome_id, experiment_id, biological_replicate, technical_replicate having c > 1 order by experiment_id';
    
    my $sth = $dbc->prepare($find_non_unique_replicate_specifications);
    $sth->execute;
    my $non_unique_replicate_specifications = $sth->fetchall_hashref('experiment_id');
    
    my @experiment_ids = sort { $a <=> $b } keys %$non_unique_replicate_specifications;
    if (@experiment_ids) {
      push @error_msg,
	"The experiments with the following experiment ids have non unique replicate specifications: (" 
	. (join ', ', @experiment_ids)
	. ") Useful sql: $find_non_unique_replicate_specifications";
    }
    
    my $find_signal_control_epigenome_mismatches = 'select * from experiment s join experiment c on (s.control_id=c.experiment_id) and s.epigenome_id != c.epigenome_id';
    
    my $sth = $dbc->prepare($find_signal_control_epigenome_mismatches);
    $sth->execute;
    my $signal_control_epigenome_mismatches = $sth->fetchall_arrayref;
    if (@$signal_control_epigenome_mismatches) {
		push @error_msg, "Signal control mismatches!";
    }
    

    
    # Check that all sequence files exist
    
    my $input_subset_files_sql = 'select input_subset_id, local_url from input_subset_tracking';
    $sth = $dbc->prepare($input_subset_files_sql);
    $sth->execute;
    
    my $local_urls = $sth->fetchall_hashref('local_url');
    my @input_subset_files = keys %$local_urls;
    
#     foreach my $current_file (@input_subset_files) {
#       if (! -e $current_file) {
# 	push @error_msg, 
# 	  "The file $current_file is specified in the input_subset ("
# 	  . $local_urls->{$current_file}->{input_subset_id}
# 	  . ") table, but it doesn't exist!"
#       }
#     }
    
    # Check for orphan result_set_inputs
    
    my $count_orphans_sql = 'select count(*) as c from result_set_input left join result_set using (result_set_id) where result_set.result_set_id is null';
    $sth = $dbc->prepare($count_orphans_sql);
    $sth->execute;
    
    my $x = $sth->fetchall_arrayref;
    my $num_orphans = $x->[0]->[0];
    
    if ($num_orphans) {
        push @error_msg, "There are $num_orphans orphan entries in the result_set_input table. They must be removed or new result_sets might have arbitrary links to input_subsets. Suggestion: delete from result_set_input where result_set_id not in (select result_set_id from result_set);"; 
    }

#     my @exported_variables = qw( PERL5LIB JAVA_HOME R_LIBS CLASSPATH );
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
    
#     my $CCAT_chr_lengths_file = $ENV{CCAT_chr_lengths_file};
#     
#     if (! $CCAT_chr_lengths_file) {
#       push @error_msg, "CCAT_chr_lengths_file has not been set!";
#     } elsif (! -e $CCAT_chr_lengths_file) {
#       push @error_msg, "$CCAT_chr_lengths_file is not a file!";
#     } else {
#       use File::Basename;
#       my $dir = dirname($CCAT_chr_lengths_file);
#       if (! -w $dir) {
# 	# This is where the sorted file gets written to.
# 	push @error_msg, "Can't write to directory $dir!";
#       }
#     }

#     my $picard_output = `java picard.cmdline.PicardCommandLine`;
#     if ($picard_output !~ /^USAGE: PicardCommandLine/) {
#       push @error_msg, qq(Can't run picard! "java picard.cmdline.PicardCommandLine".\n);
#     }
    
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
    
      die(
	"\n-------------------------------------------------------------------------------\n"
	. "Prepipeline checks have failed with the following errors:\n\n"
	. $err_string
	. "\n\nPlease fix these first before continuing with the pipeline."
	. "\n-------------------------------------------------------------------------------\n"
      );
    }
    $logger->info("Pre pipeline checks have completed successfully. You can now run the rest of the pipeline.\n");
}

1;
