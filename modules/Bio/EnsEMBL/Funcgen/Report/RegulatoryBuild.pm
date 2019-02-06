package Bio::EnsEMBL::Funcgen::Report::RegulatoryBuild;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/report.html';
  return $template;
}

sub template_dir {
  my $self = shift;
  my $template_dir = $self->template_base_dir . '/regulatory_build';
  return $template_dir;
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $previous_version_funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species . '_previous_version', 'funcgen');
  
  my $regulatory_build_adaptor
    = $funcgen_adaptor
      ->get_RegulatoryBuildAdaptor;

  my $regulatory_build_statistics_adaptor 
    = $funcgen_adaptor
      ->get_RegulatoryBuildStatisticAdaptor;

  my $regulatory_build_statistics_adaptor_previous_version 
    = $previous_version_funcgen_dba
      ->get_RegulatoryBuildStatisticAdaptor;
  
  my $genome = Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "GenomeContainer" );
  my $ref_length = $genome->get_ref_length;

  return {
    dbc      => $funcgen_adaptor->dbc,
    species  => $species,
    regulatory_build_adaptor => $regulatory_build_adaptor,
    regulatory_build_statistics_adaptor => $regulatory_build_statistics_adaptor,
    regulatory_build_statistics_adaptor_previous_version => $regulatory_build_statistics_adaptor_previous_version,
    ref_length => $ref_length,
  };
}

1;
