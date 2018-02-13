package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::FindControlExperiments;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::ExecutionPlanUtils qw (
    lock_execution_plan
    resolve_nonterminal_symbols
);

use constant {

  # Used for idr
  BRANCH_SIGNALS  => 2,
  
  # Used for aligning and fastqc
  BRANCH_CONTROLS => 3,
};

sub run {

  my $self = shift;
  my $species             = $self->param_required('species');
  my $execution_plan_list = $self->param_required('execution_plan_list');
  
  $self->say_with_header("Got " . scalar @$execution_plan_list . " execution plans.", 1);
  
  my %found_control_alignments;
  my %plan_depending_on_control;
  
  foreach my $current_execution_plan (@$execution_plan_list) {

    my $current_execution_plan_expanded = resolve_nonterminal_symbols($current_execution_plan);
    lock_execution_plan($current_execution_plan_expanded);
    
    my $alignment_plans = $current_execution_plan_expanded->{alignment};
    
    ALIGNMENT_PLAN:
    foreach my $alignment_name (keys %$alignment_plans) {
    
      my $alignment_plan = $alignment_plans->{$alignment_name};
      
#       if (! exists $alignment_plan->{is_control}) {
#         die(Dumper($alignment_plan));
#       }
      
      my $want_this
        =    
             ( $alignment_plan->{name}     ne 'No control'         )
          && ( $alignment_plan->{is_control}                       )
          && ( $alignment_plan->{analysis} eq 'remove_duplicates'  )
      ;
      if (! $want_this) {
        next ALIGNMENT_PLAN;
      }
      $found_control_alignments{$alignment_name} = $alignment_plan;
      
      if (! exists $plan_depending_on_control{$alignment_name}) {
        $plan_depending_on_control{$alignment_name} = []
      }
      push $plan_depending_on_control{$alignment_name}, $current_execution_plan;
    }
  }

  my @all_control_alignments = keys %found_control_alignments;
  
  foreach my $control_alignment (@all_control_alignments) {
  
    $self->dataflow_output_id(
      {
        'execution_plan' => $found_control_alignments{$control_alignment},
        'experiment'     => $found_control_alignments{$control_alignment}->{from_experiment},
        'species'        => $species,
      }, 
      BRANCH_CONTROLS
    );
    
    my $execution_plans_waiting_for_that_control
      = $plan_depending_on_control{$control_alignment};
    
    $self->dataflow_output_id(
      {
        'execution_plan_list' => $execution_plans_waiting_for_that_control,
        'species'             => $species,
      }, 
      BRANCH_SIGNALS
    );
  }
}

1;
