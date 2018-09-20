package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedSwemblRuns;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use constant {
  BRANCH_OUTPUT => 2,
};

sub run {

  my $self           = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($execution_plan);
  lock_execution_plan($plan_expanded);
  
  my $idr_plan = $plan_expanded->{idr};
  
  my $strategy = $idr_plan->{strategy};
  if ($strategy eq SKIP_IDR) {
    return;
  }

  my $alignment_replicates = $idr_plan->{alignment_replicates};
  
  foreach my $alignment_replicate (@$alignment_replicates) {
    $self->dataflow_output_id( 
      {
        'execution_plan' => $alignment_replicate,
        'species'        => $species,
      }, 
      BRANCH_OUTPUT
    );
  }
}

1;
