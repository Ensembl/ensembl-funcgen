package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::NoIdrComputePeakCallingInputs;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_ALIGN_ALL  => 2,
  BRANCH_CALL_PEAKS => 4,
};

sub run {

  my $self = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');
  
  my $alignment_plan = $execution_plan
    ->{call_peaks}
    ->{alignment}
    ->{source}
    ;

  my $call_peaks = $execution_plan
    ->{call_peaks};

  $self->dataflow_output_id(
    {
      'plan'    => $alignment_plan,
      'species' => $species,
    }, 
    BRANCH_ALIGN_ALL
  );
  $self->dataflow_output_id(
    {
      'plan'    => $call_peaks,
      'species' => $species,
    }, 
    BRANCH_CALL_PEAKS
  );
}

1;
