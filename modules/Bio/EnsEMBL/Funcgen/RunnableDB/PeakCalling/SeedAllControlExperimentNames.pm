package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllControlExperimentNames;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use constant {
  BRANCH_OUTPUT => 2,
};

sub run {

  my $self    = shift;
  my $species = $self->param_required('species');
  
  my $execution_plan_adaptor = Bio::EnsEMBL::Registry
      ->get_adaptor(
          $species, 
          'funcgen', 
          'ExecutionPlan'
      );

  my $execution_plan_list = [ map { $_->execution_plan_deserialised } @{$execution_plan_adaptor->fetch_all} ];
  
  my %unique_control_experiments_hash;
  foreach my $execution_plan (@$execution_plan_list) {
    my $control_experiment_name = $execution_plan->{meta_data}->{control_experiment};
    $unique_control_experiments_hash{$control_experiment_name} = 1;
  }

  my @control_experiment_names = keys %unique_control_experiments_hash;

  CONTROL_EXPERIMENT_NAME:
  foreach my $control_experiment_name (@control_experiment_names) {
  
    if ($control_experiment_name eq NA) {
      next CONTROL_EXPERIMENT_NAME;
    }
  
    $self->dataflow_output_id( 
      {
        'experiment' => $control_experiment_name,
        'species'    => $species,
      }, 
      BRANCH_OUTPUT
    );
  }
}

1;
