package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::PlanMetadataBuilder;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub create_execution_plan_meta_data {

  my $self = shift;
  my $species = shift;
  my $experiment = shift;
  
  my $feature_type = $experiment->feature_type;
  my $feature_type_name = $feature_type->name;
  
  my $feature_type_creates_broad_peaks_legible;
  
  if ($feature_type->_creates_broad_peaks) {
    $feature_type_creates_broad_peaks_legible = 'yes'
  } else {
    $feature_type_creates_broad_peaks_legible = 'no'
  }
  
  my $meta_data = {
    experiment                       => $experiment->name,
    feature_type                     => $feature_type_name,
    feature_type_creates_broad_peaks => $feature_type_creates_broad_peaks_legible,
    replicate_configurations         => $self->summarise_replicate_configurations(
      $species, 
      $experiment
    ),
  };
  return $meta_data;
}

sub summarise_replicate_configurations {

  my $self = shift;
  my $species = shift;
  my $experiment = shift;
  
  my $read_file_experimental_configuration_adaptor 
  = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ReadFileExperimentalConfiguration'
  );

  my $read_file_experimental_configuration_list 
    = $read_file_experimental_configuration_adaptor
      ->fetch_all_by_Experiment($experiment);
  
  my @summaries;
  foreach my $experimental_configuration (@$read_file_experimental_configuration_list) {
  
    my $current_summary = "("
    . $experimental_configuration->biological_replicate
    . ","
    . $experimental_configuration->technical_replicate
    . ")";
    
    push @summaries, $current_summary
  }
  return join ", ", @summaries;
}

1;
