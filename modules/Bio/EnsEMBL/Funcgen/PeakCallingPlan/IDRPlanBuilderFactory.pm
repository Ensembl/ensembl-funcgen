package Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilderFactory;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub make {

  my $self       = shift;
  my $experiment = shift;
  
  my $idr_strategy = $self->select_idr_strategy($experiment);
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::SkipIdr;
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::RunIdrOnBiologicalReplicates;
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::RunIdrOnTechnicalReplicates;
  
  my $idr_plan_builder;
  
  if ($idr_strategy eq SKIP_IDR) {
    $idr_plan_builder 
      = Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::SkipIdr
        ->new($experiment);
  }
  if ($idr_strategy eq RUN_IDR_ON_BIOLOGICAL_REPLICATES) {
    $idr_plan_builder 
      = Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::RunIdrOnBiologicalReplicates
        ->new($experiment);
  }
  if ($idr_strategy eq RUN_IDR_ON_TECHNICAL_REPLICATES) {
    $idr_plan_builder 
      = Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::RunIdrOnTechnicalReplicates
        ->new($experiment);
  }
  if (! defined $idr_strategy) {
    die("Unknown idr strategy $idr_strategy!");
  }
  return $idr_plan_builder;
}

sub select_idr_strategy {

  my $self = shift;
  my $experiment = shift;
  
  my $feature_type = $experiment->feature_type;
  
  if ($feature_type->_creates_broad_peaks) {
    return SKIP_IDR
  }

  my $number_of_biological_replicates = $experiment->count_biological_replicates;
  my $number_of_technical_replicates  = $experiment->count_technical_replicates;

  if ($number_of_biological_replicates > 1) {
    return RUN_IDR_ON_BIOLOGICAL_REPLICATES
  }
  if (
       ($number_of_biological_replicates == 1)
    && ($number_of_technical_replicates   > 1)
  ) {
      return RUN_IDR_ON_TECHNICAL_REPLICATES;
  }
  if (
       ($number_of_biological_replicates == 1)
    && ($number_of_technical_replicates  == 1)
  ) {
      return SKIP_IDR;
  }
  my $feature_type_name = $feature_type->name;
  die "Unforseen case! Feature_type: $feature_type_name for experiment: " . $experiment->name;
}

1;
