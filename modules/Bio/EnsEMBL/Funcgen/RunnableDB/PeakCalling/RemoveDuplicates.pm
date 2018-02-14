package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RemoveDuplicates;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
);

sub run {

  my $self = shift;
  
  my $data_root_dir = $self->param_required('data_root_dir');
  my $plan          = $self->param_required('plan');
  
  lock_execution_plan($plan);

  my $bam_file_no_duplicates = $plan
    ->{output}
    ->{real}
  ;

  my $bam_file_with_duplicates = $plan
    ->{input}
    ->{output}
    ->{real}
  ;
  
  eval {
    use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;
    remove_duplicates_from_bam({
      input_bam  => $data_root_dir . '/' . $bam_file_with_duplicates,
      output_bam => $data_root_dir . '/' . $bam_file_no_duplicates, 
      debug      => $self->debug,
    });
  };
  if ($@) {
    $self->throw($@);
  }

  my $cmd=qq(java picard.cmdline.PicardCommandLine ValidateSamFile INPUT=$bam_file_no_duplicates);
  $self->run_system_command($cmd);

  # Give file system time to sync
  sleep(20);
  return;
}

1;
