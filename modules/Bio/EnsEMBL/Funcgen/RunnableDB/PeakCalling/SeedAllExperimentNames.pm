package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllExperimentNames;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_OUTPUT => 2,
};

sub run {

  my $self                = shift;
  my $species             = $self->param_required('species');
  my $execution_plan_list = $self->param_required('execution_plan_list');
  
  foreach my $execution_plan (@$execution_plan_list) {
    $self->dataflow_output_id( 
      {
        'experiment' => $execution_plan->{meta_data}->{experiment},
        'species'    => $species,
      }, 
      BRANCH_OUTPUT
    );
  }
}

1;
