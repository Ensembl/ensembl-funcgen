package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::FindControlExperiments;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_IDR            => 2,
  BRANCH_ALIGN_CONTROLS => 3,
};

sub run {

  my $self = shift;
  my $species             = $self->param_required('species');
  my $execution_plan_list = $self->param_required('execution_plan_list');
  my $in_test_mode        = $self->param('in_test_mode');
  
  if ($in_test_mode) {
    $execution_plan_list = $self->reduce_execution_plan($execution_plan_list);
  }
  
  $self->say_with_header("Got " . scalar @$execution_plan_list . " execution plans.");
  
  my %seen_control_experiments;
  my %control_alignment_plan;
  
  EXECUTION_PLAN:
  foreach my $current_execution_plan (@$execution_plan_list) {
  
    my $control_experiment = $current_execution_plan
      ->{call_peaks}
      ->{control_alignment}
      ->{name}
    ;
    
    $self->say_with_header("$control_experiment");

    next EXECUTION_PLAN unless (defined $control_experiment);
    
    if (! exists $seen_control_experiments{$control_experiment}) {
      $seen_control_experiments{$control_experiment} = [];
      
      # Use first, all plans to align the control will be the same.
      #
      $control_alignment_plan{$control_experiment}
        = $current_execution_plan
          ->{call_peaks}
          ->{control_alignment}
          ->{source}
      ;
    }
    push 
      @{$seen_control_experiments{$control_experiment}},
      $current_execution_plan;
  }
  my @all_control_experiments = keys %seen_control_experiments;
  
  foreach my $control_experiment (@all_control_experiments) {
  
    $self->dataflow_output_id(
      {
        'plan'    => $control_alignment_plan{$control_experiment},
        'species' => $species,
      }, 
      BRANCH_ALIGN_CONTROLS
    );
    
    my $execution_plans_waiting_for_that_control
      = $seen_control_experiments{$control_experiment};
    
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
  
    my $current_idr_type = $execution_plan
      ->{call_peaks}
      ->{run_idr}
      ->{type}
    ;
    if (! exists $seen_idr_type{$current_idr_type}) {
    
      push @keep_execution_plan, $execution_plan;
      $seen_idr_type{$current_idr_type} = 1;
    
    }
  }
  return \@keep_execution_plan
}

1;
