package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RemoveDuplicates;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $bam_file_with_duplicates      = $self->param_required('merged_bam');
  my $full_path_to_deduplicated_bam = $self->param_required('deduplicated_bam');
  my $chunks                        = $self->param_required('chunks');
  
  use Bio::EnsEMBL::Registry;
  Bio::EnsEMBL::Registry->set_disconnect_when_inactive;
  
  eval {
    use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;
    remove_duplicates_from_bam({
      input_bam  => $bam_file_with_duplicates,
      output_bam => $full_path_to_deduplicated_bam, 
      debug      => $self->debug,
    });
  };
  if ($@) {
    $self->throw($@);
  }

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
