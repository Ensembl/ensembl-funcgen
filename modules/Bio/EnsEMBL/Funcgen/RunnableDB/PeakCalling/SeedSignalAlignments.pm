package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedSignalAlignments;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use constant {
  BRANCH_ALIGN => 2,
};

sub run {

  my $self = shift;
  my $species      = $self->param_required('species');
  my $plan         = $self->param_required('execution_plan');
  my $tempdir      = $self->param_required('tempdir_peak_calling');
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($plan);
  lock_execution_plan($plan_expanded);
  lock_execution_plan($plan);
  
  my $alignment_plans = $plan->{alignment};
  
  print Dumper($plan);
  
  my @alignment_names = keys %$alignment_plans;
  
  my @want_alignment;
  
  ALIGNMENT_PLAN:
  foreach my $alignment_name (@alignment_names) {
    my $alignment_plan = $alignment_plans->{$alignment_name};
    
    if ($alignment_plan->{name} eq NO_CONTROL_FLAG) {
      next ALIGNMENT_PLAN;
    }
    
    my $want_this_alignment = 
         $alignment_plan->{is_control} eq FALSE
      && $alignment_plan->{analysis}   eq REMOVE_DUPLICATES_ANALYSIS;
   
   next ALIGNMENT_PLAN if (! $want_this_alignment);
   
   push @want_alignment, $alignment_name;
  }

  my @alignment_plans_to_run;

  foreach my $alignment_name (@want_alignment) {  
    my $alignment = $plan_expanded->{alignment}->{$alignment_name};
    push @alignment_plans_to_run, $alignment;
  }

  foreach my $alignment_plan_to_run (@alignment_plans_to_run) {
    $self->dataflow_output_id( 
        {
          'species'        => $species,
          'execution_plan' => $alignment_plan_to_run,
        }, 
        BRANCH_ALIGN
    );
  }
  return;
}

1;
