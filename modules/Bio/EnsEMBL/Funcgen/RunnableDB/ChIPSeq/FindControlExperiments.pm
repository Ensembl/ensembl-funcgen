package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::FindControlExperiments;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::ExecutionPlanUtils qw (
    lock_execution_plan
    resolve_nonterminal_symbols
);

use constant {
  BRANCH_IDR            => 2,
  BRANCH_ALIGN_CONTROLS => 3,
};

sub run {

  my $self = shift;
  my $species             = $self->param_required('species');
  my $execution_plan_list = $self->param_required('execution_plan_list');
  my $in_test_mode        = $self->param('in_test_mode');
  $in_test_mode = 1;
  if ($in_test_mode) {
    #$execution_plan_list = $self->reduce_execution_plan($execution_plan_list);
  }
  
  $self->say_with_header("Got " . scalar @$execution_plan_list . " execution plans.", 1);
  
  my %found_control_experiments;
  my %plan_depending_on_control;
  
  foreach my $current_execution_plan (@$execution_plan_list) {

    my $current_execution_plan_expanded = resolve_nonterminal_symbols($current_execution_plan);
    lock_execution_plan($current_execution_plan_expanded);
    
    my $alignment_plans = $current_execution_plan_expanded->{alignment};
    foreach my $alignment_name (keys %$alignment_plans) {
    
      my $alignment_plan = $alignment_plans->{$alignment_name};
      
      my $want_this
        =    ( $alignment_plan->{is_control}                       )
          && ( $alignment_plan->{analysis} eq 'remove_duplicates'  )
      ;
      if ($want_this) {
        $found_control_experiments{$alignment_name} = $alignment_plan;
        
        if (! exists $plan_depending_on_control{$alignment_name}) {
          $plan_depending_on_control{$alignment_name} = []
        }
        push $plan_depending_on_control{$alignment_name}, $current_execution_plan;
      }
    }
  }

#   EXECUTION_PLAN:
#   foreach my $current_execution_plan (@$execution_plan_list) {
# 
#     use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::ExecutionPlanUtils qw (
#         lock_execution_plan
#         resolve_nonterminal_symbols
#     );
# 
#     my $current_execution_plan_expanded = resolve_nonterminal_symbols($current_execution_plan);
#     lock_execution_plan($current_execution_plan_expanded);
#     lock_execution_plan($current_execution_plan);
# 
#     my $control_experiment = $current_execution_plan_expanded
#       ->{call_peaks}
#       ->{input}
#       ->{control}
#       ->{name}
#     ;
#     
#     $self->say_with_header($control_experiment);
# 
#     next EXECUTION_PLAN unless (defined $control_experiment);
#     
#     if (! exists $found_control_experiments{$control_experiment}) {
#       $found_control_experiments{$control_experiment} = [];
#       
#       # Use first, all plans to align the control will be the same.
#       #
#       
#       my $control_alignment_plan 
#         = $current_execution_plan_expanded
#           ->{call_peaks}
#           ->{input}
#           ->{control}
#       
#       # 
#       if ($control_alignment_plan->{task} eq 'convert bam to bed') {
#         $control_alignment_plan = $control_alignment_plan->{input}
#       }
#       
#       $control_alignment_plan{$control_experiment}
#         = $control_alignment_plan
#       ;
#     }
#     push 
#       @{$found_control_experiments{$control_experiment}},
#       $current_execution_plan;
#   }
  my @all_control_experiments = keys %found_control_experiments;
  
  #die (Dumper(\%found_control_experiments));
  #die (Dumper(\%plan_depending_on_control));
  
  
  foreach my $control_experiment (@all_control_experiments) {
  
    $self->dataflow_output_id(
      {
        'execution_plan' => $found_control_experiments{$control_experiment},
        'species'        => $species,
      }, 
      BRANCH_ALIGN_CONTROLS
    );
    
    my $execution_plans_waiting_for_that_control
      = $plan_depending_on_control{$control_experiment};
    
    $self->dataflow_output_id(
      {
        'execution_plan_list' => $execution_plans_waiting_for_that_control,
        'species'             => $species,
      }, 
      BRANCH_IDR
    );
  }
}

sub reduce_execution_plan {

  my $self = shift;
  my $execution_plan_list = shift;
  
  my @keep_execution_plan;
  my %seen_idr_type;
  
  foreach my $execution_plan (@$execution_plan_list) {
    
    use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
    );

    my $execution_plan_expanded = resolve_nonterminal_symbols($execution_plan);
    lock_execution_plan($execution_plan_expanded);
    
    my $current_idr_type = $execution_plan_expanded
      ->{call_peaks}
      ->{run_idr}
      ->{strategy}
    ;
    
    if (! exists $seen_idr_type{$current_idr_type}) {
    
      push @keep_execution_plan, $execution_plan;
      $seen_idr_type{$current_idr_type} = 1;
    
    }
  }
  return \@keep_execution_plan
}

1;
