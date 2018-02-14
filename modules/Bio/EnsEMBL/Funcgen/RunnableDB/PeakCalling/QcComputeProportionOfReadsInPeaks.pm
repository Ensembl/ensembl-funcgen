package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcComputeProportionOfReadsInPeaks;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';

sub run {
  my $self = shift;

  my $peak_file     = $self->param_required('peak_file');
  my $temp_dir      = $self->param_required('frip_tempdir');
  my $peak_calling  = $self->param_required('peak_calling');
  my $bam_file      = $self->param_required('bam_file');
  my $species       = $self->param_required('species');
  
  use Bio::EnsEMBL::Utils::Logger;
  my $logger = Bio::EnsEMBL::Utils::Logger->new();
  $logger->init_log;

  if (! -e $peak_file) {
    $logger->error("The peak file $peak_file does not exist!\n");
  }
  if (! -e $bam_file) {
    $logger->error("The signal bam file $bam_file does not exist!\n");
  }

  if ($temp_dir) {
    if (-d $temp_dir) {
      $logger->info("Temporary directory $temp_dir exists.\n");
    } else {
      $logger->info("Temporary directory $temp_dir does not exist, so will create it.\n");
      make_path($temp_dir);
    }
  }

  use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd run_system_cmd );

  my $READS_IN_PEAKS_BAM_FILE = "${bam_file}.peaks.bam";

  my $current_cmd = "bedtools intersect -abam $bam_file -b $peak_file > $READS_IN_PEAKS_BAM_FILE";
  $logger->info("Running $current_cmd\n", 0, 1);
  run_system_cmd($current_cmd, 0, 1);

  $current_cmd = "samtools view -c $READS_IN_PEAKS_BAM_FILE";
  $logger->info("Running $current_cmd\n", 0, 1);
  my $num_reads_in_peaks = run_backtick_cmd($current_cmd);

  $current_cmd = "samtools view -c $bam_file";
  $logger->info("Running $current_cmd\n", 0, 1);
  my $num_reads_in_total = run_backtick_cmd($current_cmd);

  unlink($READS_IN_PEAKS_BAM_FILE);

  my $proportion_of_reads_in_peaks = $num_reads_in_peaks / $num_reads_in_total;

  $logger->info("Proportion_of_reads_in_peaks = $num_reads_in_peaks / $num_reads_in_total = $proportion_of_reads_in_peaks\n");

  my $input_id = {
  
    peak_file          => $peak_file,
    frip_tempdir       => $temp_dir,
    bam_file           => $bam_file,
    peak_calling       => $peak_calling,
    species            => $species,
    
    num_reads_in_peaks           => $num_reads_in_peaks,
    num_reads_in_total           => $num_reads_in_total,
    proportion_of_reads_in_peaks => $proportion_of_reads_in_peaks,

  };
  $self->dataflow_output_id($input_id, 2);
  return;
}

1;
