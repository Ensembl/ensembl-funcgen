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
  my $chunks        = $self->param_required('chunks');
  
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
  
  my $full_path_to_deduplicated_bam = $data_root_dir . '/' . $bam_file_no_duplicates;
  
  eval {
    use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;
    remove_duplicates_from_bam({
      input_bam  => $data_root_dir . '/' . $bam_file_with_duplicates,
      output_bam => $full_path_to_deduplicated_bam, 
      debug      => $self->debug,
    });
  };
  if ($@) {
    $self->throw($@);
  }

#   my $cmd=qq(java picard.cmdline.PicardCommandLine ValidateSamFile INPUT=$bam_file_no_duplicates);
#   $self->run_system_command($cmd);

  # Give file system time to sync
  sleep(20);

  my $cmd = "check_bam_file_has_EOF_marker.pl --bam_file $full_path_to_deduplicated_bam";
  my $has_failed = $self->run_system_command($cmd);
  if ($has_failed) {
    $self->throw("End of file marker check failed:\n" . $cmd)
  }

  my $cmd = qq(samtools index $full_path_to_deduplicated_bam);
  $has_failed = $self->run_system_command($cmd);
  if ($has_failed) {
    $self->throw("Can't index $full_path_to_deduplicated_bam")
  }
  
  my $expected_index_file = "${full_path_to_deduplicated_bam}.bai";
  if (! -e $full_path_to_deduplicated_bam) {
    $self->throw("Can't find index file ${full_path_to_deduplicated_bam}!")
  }

  foreach my $current_chunk (@$chunks) {
    
      use File::Basename qw( dirname );
      my $temp_dir = dirname($current_chunk);
      
      $self->say_with_header("Deleting $temp_dir", 1);
      
      use File::Path qw( rmtree );
      rmtree($temp_dir);
  }
  return;
}

1;
