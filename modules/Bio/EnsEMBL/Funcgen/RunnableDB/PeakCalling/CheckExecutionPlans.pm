package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CheckExecutionPlans;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
);
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use constant {
  BRANCH_WRITE_AS_ONE_LIST => 2,
};

sub run {
    my $self = shift;
    
    my $species             = $self->param_required('species');
    
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

    my @error_messages = (
        $self->check_once_a_control_always_a_control                   ($execution_plan_list),
        $self->check_alignments_with_same_name_link_to_same_experiment ($execution_plan_list),
        $self->check_read_files_are_used_only_once_in_idr              ($execution_plan_list)
     );
    
    if (@error_messages) {
        my $error_message
            = 
                "The following errors were detected:\n"
                . join ("\n", map { "  - " . $_ } @error_messages);
        $self->throw($error_message);
    }
}

sub check_once_a_control_always_a_control {

    my $self = shift;
    my $all_execution_plans = shift;

    my @error_messages;
    my %experiments_affected;
    
    my %duplicates;
    foreach my $execution_plan (@$all_execution_plans) {
        
        my $experiment_name = $execution_plan->{meta_data}->{experiment};
        
        my $alignments = $execution_plan->{alignment};
        
        my @alignment_names = keys %$alignments;
        
        ALIGNMENT:
        foreach my $alignment_name (@alignment_names) {
        
            my $alignment_plan = $alignments->{$alignment_name};
            if (! exists $alignment_plan->{analysis}) {
                print Dumper($alignment_plan);
            }
            #print Dumper($alignment_plan);
            if (! exists $alignment_plan->{analysis} || $alignment_plan->{analysis} ne REMOVE_DUPLICATES_ANALYSIS) {
                next ALIGNMENT;
            }
            
            my $is_control       = $alignment_plan->{is_control};
            my $input_is_control = $alignment_plan->{input}->{is_control};
            
            if ($is_control ne $input_is_control) {
                push 
                    @error_messages, 
                    "Inconsistent case for is_control in experiment $experiment_name:\n"
                    . Dumper($alignment_plan)
            }
        }
    }
    return @error_messages;
}

sub check_read_files_are_used_only_once_in_idr {

    my $self = shift;
    my $all_execution_plans = shift;

    my @error_messages;
    my %experiments_affected;
    
    my %duplicates;
    ALIGNMENT:
    foreach my $execution_plan (@$all_execution_plans) {
        
        my $experiment_name = $execution_plan->{meta_data}->{experiment};
        
        if (! exists $execution_plan->{idr}->{alignment_replicates}) {
            next ALIGNMENT;
        }
        
        my $alignments = $execution_plan->{idr}->{alignment_replicates};
        my %seen_read_file_names;
        
        foreach my $alignment_plan (@$alignments) {
        
            my $read_files = $alignment_plan->{input}->{input}->{read_files};
            foreach my $read_file (@$read_files) {
            
                my $type = $read_file->{type};
                
                my @names;
                
                if ($type eq SINGLE_END) {
                    push @names, $read_file->{name};
                }
                if ($type eq PAIRED_END) {
                    push @names, $read_file->{1};
                    push @names, $read_file->{2};
                }
                
                foreach my $name (@names) {
                    if (exists $seen_read_file_names{$name}) {
                        $duplicates{$name} = 1;
                        
                        if (! exists $experiments_affected{$experiment_name}) {
                            $experiments_affected{$experiment_name} = [];
                        }
                        push 
                            @{$experiments_affected{$experiment_name}},
                            $name,
                    }
                    $seen_read_file_names{$name} = 1;
                }
            }
        }
    }
    my @duplicate_read_files = keys %duplicates;
    
    foreach my $duplicate (@duplicate_read_files) {
        push @error_messages, "$duplicate is used more than once!"
    }
    foreach my $experiment (keys %experiments_affected) {
        push 
            @error_messages, 
            "The experiment $experiment has duplicate read files:\n"
            . Dumper($experiments_affected{$experiment})
    }
    
    return @error_messages;
}

sub check_alignments_with_same_name_link_to_same_experiment {

    my $self = shift;
    my $all_execution_plans = shift;

    my @error_messages;
    my %alignment_name_to_experiment_name;
    foreach my $execution_plan (@$all_execution_plans) {

        my $alignments = $execution_plan->{alignment};
        my @alignment_names = keys %$alignments;
        
        ALIGNMENT:
        foreach my $alignment_name (@alignment_names) {
            
            if (! exists $alignments->{$alignment_name}->{from_experiment}) {
                next ALIGNMENT;
            }
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
