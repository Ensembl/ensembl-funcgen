package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::MergeBamFiles;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $chunks        = $self->param_required('chunks');
  my $data_root_dir = $self->param_required('data_root_dir');
  my $species       = $self->param_required('species');
  my $merged_bam    = $self->param_required('merged_bam');
  
  use Bio::EnsEMBL::Registry;
  Bio::EnsEMBL::Registry->set_disconnect_when_inactive;
  
  use File::Basename qw( dirname basename );
  my $full_path = dirname($merged_bam);

  use File::Path qw(make_path remove_tree);
  
  $self->say_with_header("Creating directory $full_path", 1);
  make_path($full_path);

  my $full_path_to_merged_bam = $merged_bam;
  
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
