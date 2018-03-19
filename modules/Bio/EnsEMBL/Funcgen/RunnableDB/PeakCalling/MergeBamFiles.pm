package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::MergeBamFiles;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $chunks        = $self->param_required('chunks');
  my $data_root_dir = $self->param_required('data_root_dir');
  my $species       = $self->param_required('species');
  my $plan          = $self->param_required('plan');

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
      lock_execution_plan
      summarise
  );
  lock_execution_plan($plan);
  print summarise($plan);
  
  my $align_plan = $plan
    ->{input}
  ;
  my $assembly       = $align_plan->{to_assembly};
  my $bam_file_real  = $align_plan->{output}->{real};

  use File::Basename qw( dirname basename );
  my $dirname = dirname($bam_file_real);
  my $full_path = $data_root_dir . '/' . $dirname;

  use File::Path qw(make_path remove_tree);
  make_path($full_path);
  
  my $full_path_to_merged_bam = $data_root_dir . '/' . $bam_file_real;
  
  use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;# merge_bams 
  merge_bams({
    input_bams => $chunks,
    output_bam => $full_path_to_merged_bam,
    debug      => $self->debug,
  });
  
  my $cmd = "check_bam_file_has_EOF_marker.pl --bam_file $full_path_to_merged_bam";
  my $has_failed = $self->run_system_command($cmd);
  if ($has_failed) {
    $self->throw("End of file marker check failed:\n" . $cmd)
  }
  
#   my $cmd = qq(java picard.cmdline.PicardCommandLine ValidateSamFile INPUT=$full_path_to_merged_bam IGNORE_WARNINGS=true);
#   $has_failed = $self->run_system_command($cmd);
#   if ($has_failed) {
#     $self->throw("End of file marker check failed:\n" . $cmd)
#   }

  my $cmd = qq(samtools index $full_path_to_merged_bam);
  $has_failed = $self->run_system_command($cmd);
  if ($has_failed) {
    $self->throw("Can't index $full_path_to_merged_bam")
  }
  
  my $expected_index_file = "${full_path_to_merged_bam}.bai";
  if (! -e $expected_index_file) {
    $self->throw("Can't find index file ${expected_index_file}!")
  }

  # Give file system time to sync
  sleep(20);
  return;
}

1;
