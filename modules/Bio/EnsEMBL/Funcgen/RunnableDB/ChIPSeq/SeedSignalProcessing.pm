package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SeedSignalProcessing;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy;

use constant {
  BRANCH_IDR         => 2,
  BRANCH_CALL_PEAKS  => 3,
};

sub run {

  my $self                = shift;
  my $species             = $self->param_required('species');
  my $execution_plan_list = $self->param_required('execution_plan_list');
  
  foreach my $execution_plan (@$execution_plan_list) {
    $self->dataflow_output_id( 
      {
        'execution_plan' => $execution_plan,
        'species'        => $species,
      }, 
      BRANCH_IDR
    );
    $self->dataflow_output_id( 
      {
        'execution_plan' => $execution_plan,
        'species'        => $species,
      }, 
      BRANCH_CALL_PEAKS
    );
  }
}

1;
