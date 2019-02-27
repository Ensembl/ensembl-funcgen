=pod 
=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::SegmentationJobFactory

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::SegmentationJobFactory;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';

use Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator;
use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( create_species_assembly_path );

sub run {
  my $self = shift;
  
  my $tempdir                 = $self->param_required('tempdir');
  my $species                 = $self->param_required('species');
  my $data_root_dir           = $self->param_required('data_root_dir');
  my $reference_data_root_dir = $self->param_required('reference_data_root_dir');
  
  my $species_assembly_path = create_species_assembly_path($species);
  
  my $data_root_dir_species_assembly = $data_root_dir . '/' . $species_assembly_path;
  
  my $binarized_bam_dir     = $tempdir . '/binarization/';
  my $learnmodel_output_dir = $tempdir . '/learn_model/';
  
  my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'coordsystem');
  my $default_chromosome_coordsystem = $coordsystem_adaptor->fetch_by_name('chromosome');
  my $default_assembly = $default_chromosome_coordsystem->version;

  my $chromosome_length_file = 
    $reference_data_root_dir
    . '/'
    . 
      Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator
        ->new
        ->locate({
          species   => $species,
          file_type => 'chromosome_lengths_by_species_assembly',
          
          # Use male to make sure the Y chromosome is in the chomosome length
          # file.
          epigenome_gender => 'male',
          assembly  => $default_assembly
        });
  
  use Data::Dumper;
  
  my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
  my $core_dbc = $core_adaptor->dbc;

  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  my $funcgen_dbc = $funcgen_adaptor->dbc;
  
  my $segmentation_parameter_file = $tempdir . '/segmentation_parameter_file.txt';
  
  my $core_dbc_password    = $core_dbc->password;
  my $funcgen_dbc_password = $funcgen_dbc->password;
  
  # If there is no password, set it to "". These parameters areare used on 
  # the command line.
  #
  if (! defined $core_dbc_password || $core_dbc_password eq '') {
    $core_dbc_password = '""';
  }
  if (! defined $funcgen_dbc_password || $funcgen_dbc_password eq '') {
    $funcgen_dbc_password = '""';
  }

  my $input_id = {
    chromosome_length_file         => $chromosome_length_file,
    data_root_dir_species_assembly => $data_root_dir_species_assembly,
    binarized_bam_dir              => $binarized_bam_dir,
    learnmodel_output_dir          => $learnmodel_output_dir,
    assembly                       => $default_assembly,
    species                        => $species,
    
    segmentation_parameter_file => $segmentation_parameter_file,
    
    core_host     => $core_dbc->host,
    core_port     => $core_dbc->port,
    core_username => $core_dbc->username,
    core_password => $core_dbc_password,
    core_dbname   => $core_dbc->dbname,

    funcgen_host     => $funcgen_dbc->host,
    funcgen_port     => $funcgen_dbc->port,
    funcgen_username => $funcgen_dbc->username,
    funcgen_password => $funcgen_dbc_password,
    funcgen_dbname   => $funcgen_dbc->dbname,

  };
  
  $self->dataflow_output_id($input_id, 2);
  return;
}

1;
