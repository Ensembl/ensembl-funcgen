=pod 
=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::RegulatoryBuild::StableIdMappingJobFactory

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::RegulatoryBuild::StableIdMappingJobFactory;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';

use Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator;
use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( create_species_assembly_path );

use constant {
  OUTPUT_BRANCH => 2,
};


sub run {
  my $self = shift;
  
  my $species = $self->param_required('species');
  my $tempdir = $self->param_required('tempdir');
  
  my $binarized_bam_dir     = $tempdir . '/binarization/';
  my $learnmodel_output_dir = $tempdir . '/learn_model/';
  
  use Data::Dumper;
  
  my $species_suffix = '_previous_version';
  
  my $previous_version_species = $species . $species_suffix;
  
  my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
  my $core_dbc = $core_adaptor->dbc;

  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  my $funcgen_dbc = $funcgen_adaptor->dbc;

  my $previous_version_funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($previous_version_species, 'funcgen');
  my $previous_version_funcgen_dbc = $previous_version_funcgen_adaptor->dbc;

  my $core_dbc_password    = $core_dbc->password;
  my $funcgen_dbc_password = $funcgen_dbc->password;

  if (! defined $core_dbc_password || $core_dbc_password eq '') {
    $core_dbc_password = '""';
  }
  if (! defined $funcgen_dbc_password || $funcgen_dbc_password eq '') {
    $funcgen_dbc_password = '""';
  }

  my $core_db_url 
    = Bio::EnsEMBL::Hive::DBSQL::DBConnection->new(
        -dbconn => $core_dbc,
    )->url;

  my $funcgen_db_url 
    = Bio::EnsEMBL::Hive::DBSQL::DBConnection->new(
        -dbconn => $funcgen_dbc,
    )->url;
    
  my $previous_version_funcgen_db_url
    = Bio::EnsEMBL::Hive::DBSQL::DBConnection->new(
        -dbconn => $previous_version_funcgen_dbc,
    )->url;
  
  my $stable_id_prefix;

  # Using a regular expression allows this to work for things like
  # "homo_sapiens_previous_version" or "homo_sapiens_rerun"
  #
  if ($species =~ /homo_sapiens.*/ ) {
    $stable_id_prefix = 'ENSR';
  }
  if ($species eq 'mus_musculus' ) {
    $stable_id_prefix = 'ENSMUSR';
  }
  
  if (! defined $stable_id_prefix) {
    die;
  }
  
  my $regulatory_features_file = $tempdir . '/' . 'regulatory_features.bed';
  my $regulatory_features_previous_version_bed_file = $tempdir . '/' . 'regulatory_features_previous_version.bed';
  my $overlaps_bed_file        = $tempdir . '/' . 'regulatory_features_overlaps.bed';
  my $stable_id_mapping_file   = $tempdir . '/' . 'stable_id_mapping.tsv';
  my $mapping_report           = $tempdir . '/' . 'mapping_report.pl';
  
  my $input_id = {
    regulatory_features_bed_file => $regulatory_features_file,
    stable_id_prefix             => $stable_id_prefix,
    db_url_funcgen               => $funcgen_db_url,
    db_url_funcgen_old           => $previous_version_funcgen_db_url,
    regulatory_features_previous_version_bed_file => $regulatory_features_previous_version_bed_file,
    overlaps_bed_file            => $overlaps_bed_file,
    stable_id_mapping_file       => $stable_id_mapping_file,
    species                      => $species,
    mapping_report               => $mapping_report,
    
#     core_host     => $core_dbc->host,
#     core_port     => $core_dbc->port,
#     core_username => $core_dbc->username,
#     core_password => $core_dbc_password,
#     core_dbname   => $core_dbc->dbname,
# 
#     funcgen_host     => $funcgen_dbc->host,
#     funcgen_port     => $funcgen_dbc->port,
#     funcgen_username => $funcgen_dbc->username,
#     funcgen_password => $funcgen_dbc_password,
#     funcgen_dbname   => $funcgen_dbc->dbname,
  };
  
  $self->dataflow_output_id($input_id, OUTPUT_BRANCH);
  return;
}

1;
