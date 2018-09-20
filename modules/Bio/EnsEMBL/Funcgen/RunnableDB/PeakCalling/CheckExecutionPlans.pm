package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CheckExecutionPlans;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
);
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Checks qw ( :all );

sub run {
    my $self    = shift;
    my $species = $self->param_required('species');
    
    my $execution_plan_adaptor = Bio::EnsEMBL::Registry
        ->get_adaptor(
            $species, 
            'funcgen', 
            'ExecutionPlan'
        );

    my $execution_plan_list = [
    
        map { lock_execution_plan($_); $_     }
        map { $_->execution_plan_deserialised }
        
        @{ $execution_plan_adaptor->fetch_all }
    ];
    
    $self->say_with_header("Got " . scalar @$execution_plan_list . " execution plans.", 1);

    my @error_messages = run_all_checks(@$execution_plan_list);
    
    if (@error_messages) {
        my $error_message
            = 
                "The following errors were detected:\n"
                . join ("\n", map { "  - " . $_ } @error_messages);
        $self->throw($error_message);
    }
    $self->say_with_header("All tests passed successfully.", 1);
    return;
}

1;
