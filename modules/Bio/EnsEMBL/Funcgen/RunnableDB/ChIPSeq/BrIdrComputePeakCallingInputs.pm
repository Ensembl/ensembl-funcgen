package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::BrIdrComputePeakCallingInputs;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_ALIGN_ALL_READ_FILES     => 3,
  BRANCH_ALIGN_REPLICATES_FOR_IDR => 2,
  BRANCH_RUN_IDR_PEAK_CALLING     => 4,
};

sub run {

  my $self = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');
  
  my $read_files = $execution_plan
    ->{call_peaks}
    ->{alignment}
    ->{source}
  ;
  my $run_idr = $execution_plan
    ->{call_peaks}
    ->{run_idr}
  ;
  my $alignment_replicates = $run_idr
    ->{alignment_replicates}
  ;

  $self->dataflow_output_id(
    {
        'plan'    => $read_files,
        'species' => $species,
    }, 
    BRANCH_ALIGN_ALL_READ_FILES
  );
  
  foreach my $alignment_replicate (@$alignment_replicates) {
  
    $self->dataflow_output_id(
        {
            'plan' => $alignment_replicate,
            'species' => $species,
        }, 
        BRANCH_ALIGN_REPLICATES_FOR_IDR
    );
    $self->dataflow_output_id(
        {
            'plan' => $alignment_replicate,
            'species' => $species,
        }, 
        BRANCH_RUN_IDR_PEAK_CALLING
    );
  }
}

1;
