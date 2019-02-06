package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::FindControlExperiments;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
    resolve_nonterminal_symbols
);

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use constant {

  BRANCH_SIGNALS                  => 2,
  BRANCH_CONTROLS                 => 3,
  BRANCH_SIGNALS_WITHOUT_CONTROLS => 4,

};

sub run {

  my $self = shift;
  my $species             = $self->param_required('species');
  #my $execution_plan_list = $self->param_required('execution_plan_list');
  my $execution_plan_adaptor = Bio::EnsEMBL::Registry
      ->get_adaptor(
          $species, 
          'funcgen', 
          'ExecutionPlan'
      );
  my $execution_plan_list = [ map { $_->execution_plan_deserialised } @{$execution_plan_adaptor->fetch_all} ];

  $self->say_with_header("Got " . scalar @$execution_plan_list . " execution plans.", 1);
  
  my %found_control_alignments;
  my %plan_depending_on_control;
  my @signal_experiments_without_controls;
  
  EXECUTION_PLAN:
  foreach my $current_execution_plan (@$execution_plan_list) {
  
    my $meta_data = $current_execution_plan->{meta_data};
    
     if (
          ( $meta_data->{experiment_is_control}  ne TRUE  )
       && ( $meta_data->{experiment_has_control} eq FALSE )
     ) {
       push 
         @signal_experiments_without_controls, 
         $current_execution_plan;
 
       next EXECUTION_PLAN;
     }

    my $current_execution_plan_expanded;
    
    #eval {
      $current_execution_plan_expanded = resolve_nonterminal_symbols($current_execution_plan);
#     };
#     if ($@) {
#       warn("Couldn't use an execution plan!");
#     }
    
    lock_execution_plan($current_execution_plan_expanded);
    
    my $alignment_plans = $current_execution_plan_expanded->{alignment};
    
    ALIGNMENT_PLAN:
    foreach my $alignment_name (keys %$alignment_plans) {
    
      my $alignment_plan = $alignment_plans->{$alignment_name};
      
      my $want_this
        =    
             ( $alignment_plan->{name}       ne NO_CONTROL_FLAG             )
          && ( $alignment_plan->{is_control} eq TRUE                        )
          && ( $alignment_plan->{analysis}   eq REMOVE_DUPLICATES_ANALYSIS  )
      ;
      if (! $want_this) {
        next ALIGNMENT_PLAN;
      }
      $found_control_alignments{$alignment_name} = $alignment_plan;
      
      if (! exists $plan_depending_on_control{$alignment_name}) {
        $plan_depending_on_control{$alignment_name} = []
      }
      push @{$plan_depending_on_control{$alignment_name}}, $current_execution_plan;
    }
  }

  foreach my $execution_plan (@signal_experiments_without_controls) {
    $self->dataflow_output_id(
      {
        'execution_plan' => $execution_plan,
        'species'        => $species,
      }, 
      BRANCH_SIGNALS_WITHOUT_CONTROLS
    );
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
 
    for my $execution_plan (@$execution_plans_waiting_for_that_control) {
    $self->dataflow_output_id(
      {   
        'execution_plan' => $execution_plan,
        'species'             => $species,
      },  
      BRANCH_SIGNALS
    );  
    }
#     $self->dataflow_output_id(
#       {
#         'execution_plan_list' => $execution_plans_waiting_for_that_control,
#         'species'             => $species,
#       }, 
#       BRANCH_SIGNALS
#     );
  }
}

1;
