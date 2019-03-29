package Bio::EnsEMBL::Funcgen::PeakCallingPlan::Checks;

use strict;
use Data::Dumper;

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use base qw( Exporter );
use vars qw( @EXPORT_OK );
use strict;

our @EXPORT_OK = qw(

  run_all_checks

  check_once_a_control_always_a_control
  check_read_files_are_used_only_once_in_idr
  check_alignments_with_same_name_link_to_same_experiment

);

our %EXPORT_TAGS = (all => \@EXPORT_OK);

sub run_all_checks {

  my @execution_plan_list = @_;

  my @error_messages;

  push @error_messages, map { check_once_a_control_always_a_control      ($_) } @execution_plan_list;
  push @error_messages, map { check_read_files_are_used_only_once_in_idr ($_) } @execution_plan_list;
  push @error_messages, check_alignments_with_same_name_link_to_same_experiment ( \@execution_plan_list );

  return @error_messages;
}

sub check_once_a_control_always_a_control {

    my $execution_plan = shift;

    my @error_messages;

    my $experiment_name = $execution_plan->{meta_data}->{experiment};
    my $alignments      = $execution_plan->{alignment};
    my @alignment_names = keys %$alignments;
    
    ALIGNMENT:
    foreach my $alignment_name (@alignment_names) {
    
        my $alignment_plan = $alignments->{$alignment_name};
        
        if (
            ! exists $alignment_plan->{analysis} 
            || $alignment_plan->{analysis} ne REMOVE_DUPLICATES_ANALYSIS
            ) {
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
    return @error_messages;
}

sub check_read_files_are_used_only_once_in_idr {

    my $execution_plan = shift;

    my @error_messages;
    my %experiments_affected;
    
    my %duplicates;
    my $experiment_name = $execution_plan->{meta_data}->{experiment};
    
    if (! exists $execution_plan->{idr}->{alignment_replicates}) {
        return @error_messages;
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
