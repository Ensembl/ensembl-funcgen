package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedBigWigControlJobs;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use constant {
  BRANCH_TO_BIGWIG => 2,
};

sub run {

  my $self    = shift;
  my $species = $self->param_required('species');
  my $plan    = $self->param_required('execution_plan');
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($plan);
  lock_execution_plan($plan_expanded);
  lock_execution_plan($plan);
  
  print Dumper($plan_expanded);
  
  my $signal_plans = $plan->{signal};
  
  my @signal_names = keys %$signal_plans;
  
  my @want_signal;
  
  SIGNAL_PLAN:
  foreach my $name (@signal_names) {
    my $signal_plan = $signal_plans->{$name};
    
    if ($signal_plan->{name} eq NO_CONTROL_FLAG) {
      next SIGNAL_PLAN;
    }

    my $want_this = 
         $signal_plan->{is_control} eq TRUE
      && $signal_plan->{analysis}   eq CONVERT_BAM_TO_BIGWIG_ANALYSIS;
   
   next SIGNAL_PLAN if (! $want_this);
   
   push @want_signal, $name;
  }

  my @signal_plans_to_run;

  foreach my $name (@want_signal) {  
    my $signal_plan = $plan_expanded->{signal}->{$name};
    push @signal_plans_to_run, $signal_plan;
  }

  foreach my $plan_to_run (@signal_plans_to_run) {
    $self->dataflow_output_id( 
        {
          'species'        => $species,
          'execution_plan' => $plan_to_run,
        }, 
        BRANCH_TO_BIGWIG
    );
  }
  return;
}

1;
