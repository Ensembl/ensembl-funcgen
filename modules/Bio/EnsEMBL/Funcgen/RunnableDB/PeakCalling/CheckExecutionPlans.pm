package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CheckExecutionPlans;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
    resolve_nonterminal_symbols
);

use constant {
  BRANCH_WRITE_AS_ONE_LIST         => 2,
#   BRANCH_WRITE_EACH_EXECUTION_PLAN => 3,
};

sub run {
    my $self = shift;
    
    my $species             = $self->param_required('species');
    my $execution_plan_list = $self->param_required('execution_plan_list');

    $self->say_with_header("Got " . scalar @$execution_plan_list . " execution plans.", 1);

    my @error_messages = $self->check_alignments_with_same_name_link_to_same_experiment($execution_plan_list);
    
    if (@error_messages) {
        my $error_message
            = 
                "The following errors were detected:\n"
                . join ("\n", map { "  - " . $_ } @error_messages);
        $self->throw($error_message);
    }

    $self->dataflow_output_id(
        {
            'species'             => $species,
            'execution_plan_list' => $execution_plan_list,
        }, 
        BRANCH_WRITE_AS_ONE_LIST
    );
#     foreach my $execution_plan (@$execution_plan_list) {
#       $self->dataflow_output_id(
#           {
#               'species'        => $species,
#               'execution_plan' => $execution_plan,
#           }, 
#           BRANCH_WRITE_EACH_EXECUTION_PLAN
#       );
#     }
}

sub check_alignments_with_same_name_link_to_same_experiment {

    my $self = shift;
    my $all_execution_plans = shift;

    my @error_messages;
    my %alignment_name_to_experiment_name;
    foreach my $execution_plan (@$all_execution_plans) {

        my $alignments = $execution_plan->{alignment};
        my @alignment_names = keys %$alignments;
        
        foreach my $alignment_name (@alignment_names) {
        
            my $experiment_name = $alignments->{$alignment_name}->{from_experiment};
        
            if (! exists $alignment_name_to_experiment_name{$alignment_name}) {
                $alignment_name_to_experiment_name{$alignment_name} = $experiment_name;
            }
            if ($alignment_name_to_experiment_name{$alignment_name} ne $experiment_name) {
                push 
                    @error_messages, 
                    "The alignment name $alignment_name created for the experiment " 
                    . $alignment_name_to_experiment_name{$alignment_name} 
                    . " is also already being used for the experiment $experiment_name";
            }
        }
    }
    return @error_messages;
}

1;
