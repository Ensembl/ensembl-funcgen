package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::WriteBigWig;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self           = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');
  my $data_root_dir  = $self->param_required('data_root_dir');

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($execution_plan);
  #lock_execution_plan($plan_expanded);
  lock_execution_plan($execution_plan);

   print Dumper($plan_expanded);
#   die;

  my $bam_file = $execution_plan
    ->{input}
    ->{output}
    ->{real};

  my $bigwig_file = $execution_plan
    ->{output}
    ->{real};
  
  my $epigenome_gender = $execution_plan
    ->{input}
    ->{input}
    ->{to_gender};
  
  my $assembly = $execution_plan
    ->{input}
    ->{input}
    ->{to_assembly};
  
  my $bam_file_full_path    = "$data_root_dir/$bam_file";
  my $bigwig_file_full_path = "$data_root_dir/$bigwig_file";

  use Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator->new;
  
  my $chromosome_lengths_relative = $bwa_index_locator->locate({
    species          => $species,
    epigenome_gender => $epigenome_gender,
    assembly         => $assembly,
    file_type        => 'chromosome_lengths_by_species_assembly',
  });
  my $reference_data_root_dir = $self->param('reference_data_root_dir');

  my $chromosome_length_file = $reference_data_root_dir . '/' . $chromosome_lengths_relative;
  
  $self->say_with_header("bam_file    = $data_root_dir/$bam_file");
  $self->say_with_header("bigwig_file = $data_root_dir/$bigwig_file");
  
  my $error_occurred;
  my $cmd;
  
  $cmd = "samtools index $bam_file_full_path";
  $error_occurred = $self->run_system_command($cmd);
  if ($error_occurred) {
    $self->throw("The following command failed:\n$cmd");
  }

  my $cmd = "samtools idxstats $bam_file_full_path | awk '{total = total + \$3} END { print total }' >&2";
  (
    $error_occurred,
    my $total_mapped
  ) = 
    $self->run_system_command(
        $cmd,
        {
            use_bash_pipefail => 1,
        }
    );
  if ($error_occurred) {
    $self->throw("The following command failed:\n$cmd");
  }
  my $wiggle_tools_cmd = 'mean scale '.(10**9 / $total_mapped).' '.$bam_file_full_path;
  
  $self->say_with_header("total_mapped = $total_mapped");
  
  use File::Basename qw( dirname basename );
  my $dirname = dirname($bigwig_file_full_path);

  use File::Path qw(make_path remove_tree);
  make_path($dirname);
  
  my $cmd_wiggletools = "wiggletools write - ". $wiggle_tools_cmd;
  my $cmd_wigToBigWig = 'wigToBigWig -fixedSummaries stdin ' . $chromosome_length_file . ' ' . $bigwig_file_full_path;
  
  my $cmd = $cmd_wiggletools . ' | ' . $cmd_wigToBigWig;
   
  $self->say_with_header($cmd);
  $error_occurred = $self->run_system_command($cmd);
  
  # If lsf kills the command, give it time to kill the worker as well so 
  # eHive handles this properly.
  #
  sleep(20);
  if ($error_occurred) {
    $self->throw("The following command failed:\n$cmd");
  }
  if (! -e $bigwig_file_full_path) {
    $self->throw(
      "The following command failed:\n\n"
      . "$cmd\n\n"
      . "The following file should have been created:\n\n"
      . "$bigwig_file_full_path\n\n"
      . "but it doesn't exist."
    );
  }
}

1;
