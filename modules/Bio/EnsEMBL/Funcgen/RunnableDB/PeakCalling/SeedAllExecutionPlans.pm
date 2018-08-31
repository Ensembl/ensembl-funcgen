package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllExecutionPlans;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

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
  
  foreach my $execution_plan (@$execution_plan_list) {
    $self->dataflow_output_id( 
      {
        'execution_plan' => $execution_plan,
        'species'        => $species,
      }, 
      BRANCH_OUTPUT
    );
  }
}

1;
