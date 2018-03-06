package Bio::EnsEMBL::Funcgen::PeakCallingPlan::PlanMetadataBuilder;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

sub set_meta_data {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('meta_data ', undef, $obj);
}

sub get_meta_data  { 
  return shift->_generic_get('meta_data ');
}

sub construct {

  my $self  = shift;
  my $param = shift;
  my $species    = $param->{species};
  my $experiment = $param->{experiment};
  
  my $feature_type = $experiment->feature_type;
  my $feature_type_name = $feature_type->name;
  
  my $feature_type_creates_broad_peaks_legible;
  
  if ($feature_type->_creates_broad_peaks) {
    $feature_type_creates_broad_peaks_legible = TRUE
  } else {
    $feature_type_creates_broad_peaks_legible = FALSE
  }
  
  my $experiment_has_control;
  if ($experiment->get_control) {
    $experiment_has_control = TRUE;
  } else {
    $experiment_has_control = FALSE;
  }
  
  my $experiment_is_control;
  if ($experiment->is_control) {
    $experiment_is_control = TRUE;
  } else {
    $experiment_is_control = FALSE;
  }
  
  my $control_experiment_name;
  if (
         ($experiment_is_control  eq TRUE) 
      || ($experiment_has_control eq FALSE)
     ) {
    $control_experiment_name = NA;
  } else {
    my $control_experiment = $experiment->get_control;
    
    if (! defined $control_experiment) {
      die("Can't get control for experiment " . $experiment->name);
    }
    
    $control_experiment_name = $control_experiment->name;
  }
  
  my $meta_data = {
    experiment                       => $experiment->name,
    feature_type                     => $feature_type_name,
    feature_type_creates_broad_peaks => $feature_type_creates_broad_peaks_legible,
    experiment_has_control           => $experiment_has_control,
    experiment_is_control            => $experiment_is_control,
    control_experiment               => $control_experiment_name,
    replicate_configurations         => $self->summarise_replicate_configurations(
      $species, 
      $experiment
    ),
  };
  $self->set_meta_data($meta_data);
  return;
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
  EXPERIMENTAL_CONFIGURATION:
  foreach my $experimental_configuration (@$read_file_experimental_configuration_list) {
    
    # Report paired end pairs only once
    if ($experimental_configuration->paired_end_tag == 2) {
      next EXPERIMENTAL_CONFIGURATION;
    }
    
    # If there are multiples, only report the first.
    if ($experimental_configuration->multiple != 1) {
      next EXPERIMENTAL_CONFIGURATION;
    }
  
    my $current_summary = "("
    . "BR" . $experimental_configuration->biological_replicate
    . ","
    . "TR" . $experimental_configuration->technical_replicate
    . ")";
    
    push @summaries, $current_summary
  }
  return join ", ", @summaries;
}

1;
