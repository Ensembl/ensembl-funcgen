package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Director;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::PeakCallingStrategy;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRPlanBuilderFactory;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::PlanMetadataBuilder;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignAllPlanBuilder;

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);


sub construct_execution_plan {

  my $self  = shift;
  my $param = shift;
  
  my $species                = $param->{species};
  my $experiment             = $param->{experiment};
  my $directory_name_builder = $param->{directory_name_builder};
  my $assembly               = $param->{assembly};

  my $idr_plan_builder_factory
    = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRPlanBuilderFactory
      ->new;
  
  my $idr_plan_builder = $idr_plan_builder_factory
    ->make_idr_plan_builder($experiment);
    
  my $idr_plan = $idr_plan_builder
    ->construct($param);
  
  my $peak_calling_strategy = $self->select_peak_calling_strategy($experiment);
  
  my $align_all_plan_builder 
    = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignAllPlanBuilder
      ->new;
  
  my $align_all_read_files_for_experiment_plan 
    = $align_all_plan_builder
      ->construct($param);
  
  my %control_param = %$param;
  $control_param{experiment} = $param->{experiment}->get_control;
  
  my $align_all_read_files_for_control_plan 
    = $align_all_plan_builder
      ->construct(\%control_param);
  
  my $plan_meta_data_builder 
    = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::PlanMetadataBuilder
      ->new;
   
  my $meta_data = $plan_meta_data_builder
    ->create_execution_plan_meta_data(
      $species, 
      $experiment
    );
  
  my $execution_plan = {
    meta_data  => $meta_data,
    call_peaks => {
      alignment             => {
        source => $align_all_read_files_for_experiment_plan,
        name   => $align_all_read_files_for_experiment_plan
          ->{remove_duplicates}
          ->{name},
      },
      run_idr               => $idr_plan,
      control_alignment     => {
        source => $align_all_read_files_for_control_plan,
        name   => $align_all_read_files_for_control_plan
          ->{remove_duplicates}
          ->{name},
      },
      peak_calling_strategy => $peak_calling_strategy,
    }
  };
  return $execution_plan;
}

sub select_peak_calling_strategy {

  my $self = shift;
  my $experiment = shift;
  
  my $feature_type      = $experiment->feature_type;

  if ($feature_type->_creates_broad_peaks) {
    return Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::PeakCallingStrategy->CALL_BROAD_PEAKS;
  }
  return Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::PeakCallingStrategy->CALL_NARROW_PEAKS
}

1;
