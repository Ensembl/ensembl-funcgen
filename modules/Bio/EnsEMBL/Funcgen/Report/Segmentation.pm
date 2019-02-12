package Bio::EnsEMBL::Funcgen::Report::Segmentation;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/report.html';
  return $template;
}

sub template_dir {
  my $self = shift;
  my $template_dir = $self->template_base_dir . '/segmentation';
  return $template_dir;
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $segmentation_statistic_adaptor = $funcgen_adaptor->get_SegmentationStatisticAdaptor;
  my $segmentation_adaptor           = $funcgen_adaptor->get_SegmentationAdaptor;

  my @segmentation_assignments = qw(
    ctcf
    dead
    distal
    gene
    poised
    proximal
    repressed
    tss
    weak
  );

  my $segmentation_state_emission_adaptor = $funcgen_adaptor->get_SegmentationStateEmissionAdaptor;

  my $segmentation_state_emissions = $segmentation_state_emission_adaptor->fetch_all;

  my $segmentation_state_emissions_sorted = [ sort { $a->state <=> $b->state } @$segmentation_state_emissions ];
  
  return {

        segmentation_statistic_adaptor 
          => $segmentation_statistic_adaptor,

        segmentation_adaptor 
          => $segmentation_adaptor,

        segmentation_assignments
          => \@segmentation_assignments,
          
        segmentation_state_emissions => $segmentation_state_emissions_sorted,
        
        dbc => $funcgen_adaptor->dbc,
        
        emission_to_rgb => sub {
            my $number = shift;
            if ($number < 0.1) {
                return '#FFFFFF';
            }
            if ($number < 0.2) {
                return '#DDDDFF';
            }
            if ($number < 0.3) {
                return '#BBBBFF';
            }
            if ($number < 0.4) {
                return '#AAAAFF';
            }
            if ($number < 0.5) {
                return '#9999FF';
            }
            if ($number < 0.6) {
                return '#8888FF';
            }
            if ($number < 0.7) {
                return '#7777FF';
            }
            if ($number < 0.8) {
                return '#6666FF';
            }
            return '#5555FF';
        },
    };
}

1;
