package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcProportionOfReadsInPeaksJobFactory;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';

sub run {
  my $self = shift;
  my $peak_calling_name = $self->param_required('peak_calling');
  my $tempdir           = $self->param_required('tempdir');
  my $species           = $self->param_required('species');
  my $data_root_dir     = $self->param_required('data_root_dir');
  my $reg_conf          = $self->param_required('reg_conf');
  
  my $peak_calling_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'PeakCalling'
    );
  my $peak_calling = $peak_calling_adaptor->fetch_by_name($peak_calling_name);

  my $peak_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Peak'
    );
  
  my $this_runs_tempdir = $tempdir . '/' . $peak_calling_name;
  
  use File::Path qw( make_path );
  make_path($this_runs_tempdir);
 
  my $bed_file = $this_runs_tempdir . '/' . $peak_calling_name . '.bed';
  open my $bed_fh, '>', $bed_file;
  $peak_adaptor->_bulk_export_to_bed_by_PeakCalling($peak_calling, $bed_fh);
  $bed_fh->close;
  
  my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'coordsystem');

  my $default_chromosome_coordsystem = $coordsystem_adaptor->fetch_by_name('chromosome');
  my $default_assembly = $default_chromosome_coordsystem->version;
  
  my $signal_alignment = $peak_calling->fetch_signal_Alignment;
  
  my $signal_bam_data_file = $signal_alignment->fetch_bam_DataFile;
  my $signal_bam_file_name = $signal_bam_data_file->path;
  
  my $signal_bam_file  = $data_root_dir . '/' . $species . '/' . $default_assembly . '/' . $signal_bam_file_name;
  
  my $input_id = {
  
    peak_file          => $bed_file,
    frip_tempdir       => $this_runs_tempdir,
    bam_file           => $signal_bam_file,
    peak_calling       => $peak_calling_name,
    species            => $species,

  };
  $self->dataflow_output_id($input_id, 2);
  return;
}

1;
