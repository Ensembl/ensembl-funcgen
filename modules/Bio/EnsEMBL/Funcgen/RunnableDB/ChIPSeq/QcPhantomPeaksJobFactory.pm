package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::QcPhantomPeaksJobFactory;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';

sub run {
  my $self = shift;

  my $species         = $self->param_required('species');
  my $alignment_name  = $self->param_required('alignment');
  my $data_root_dir   = $self->param_required('data_root_dir');
  my $tempdir         = $self->param_required('tempdir');

  my $alignment_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Alignment'
    );
  my $alignment = $alignment_adaptor->fetch_by_name($alignment_name);
  
  my $epigenome_production_name 
    = $alignment
      ->fetch_Experiment
      ->epigenome
      ->production_name;
  
  my $bam_file = $alignment->fetch_bam_DataFile->relative_ftp_site_path;
  my $bam_file_full_path = $data_root_dir . '/' . $bam_file;
  
  if (! -e $bam_file_full_path) {
    die("Can't find bam file: $bam_file_full_path");
  }
  
  my $temp_dir_for_this_run = "$tempdir/$epigenome_production_name/$alignment_name";
  
  use File::Path qw( make_path );
  make_path($temp_dir_for_this_run);
  
  use File::Basename;
  (my $bam_file_base_name,  my $signal_bam_directory)  = fileparse($bam_file);

  my $signal_phantom_peaks_file  = "$temp_dir_for_this_run/${bam_file_base_name}.phantom_peaks.txt";
  
  my $input_id = {
      # Directory into which the bam files will be copied
      tempdir               => $temp_dir_for_this_run,
      alignment_name        => $alignment_name,
      phantom_peak_out_file => $signal_phantom_peaks_file,
      bam_file              => $bam_file_full_path,
      
      species => $species,
  };
  
  $self->dataflow_output_id($input_id, 2);
  return;
}

1;


