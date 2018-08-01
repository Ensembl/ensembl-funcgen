=pod 
=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::SeedLoadAssignmentJobs

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::SeedLoadAssignmentJobs;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';
use Hash::Util qw( lock_keys );

use constant {
  BRANCH_TO_SEED_JOBS => 2,
};

sub run {
  my $self = shift;
  
  my $segmentation_lists = $self->param_required('segmentation_lists');
  
  my @superclasses = keys %$segmentation_lists;
  
  foreach my $superclass (@superclasses) {
  
    my $class_to_segmentation_meta_data = $segmentation_lists->{$superclass};
    my @classes = keys %$class_to_segmentation_meta_data;
    
    foreach my $class (@classes) {
    
        my $segmentation_meta_data = $class_to_segmentation_meta_data->{$class};
        
        lock_keys(%$segmentation_meta_data);

        my $segmentation_name     = $segmentation_meta_data->{segmentation_name};
        my $learn_model_directory = $segmentation_meta_data->{learn_model_directory};

        my $hive_job = {
            learn_model_directory => $segmentation_meta_data->{learn_model_directory},
            segmentation_name     => $segmentation_meta_data->{segmentation_name},
            celltable_file        => $segmentation_meta_data->{celltable_file},
            class                 => $class,
            superclass            => $superclass,
        };
        use Data::Dumper;
        print Dumper($hive_job);

        $self->dataflow_output_id(
            $hive_job, 
            BRANCH_TO_SEED_JOBS
        );
    }
  }
  return;
}

1;
