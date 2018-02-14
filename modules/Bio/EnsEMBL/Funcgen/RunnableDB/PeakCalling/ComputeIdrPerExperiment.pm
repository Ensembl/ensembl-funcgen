package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::ComputeIdrPerExperiment;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_OUTPUT => 2,
};

sub run {

  my $self = shift;
  
  my $execution_plan       = $self->param_required('execution_plan');
  my $species              = $self->param_required('species');
  
  # Comes from accu, so autoflow wouldn't work here.
  #
  my $permissive_peak_calling = $self->param_required('permissive_peak_calling');
  
  print Dumper($permissive_peak_calling);
  
#   use Hash::Util qw( lock_keys );
#   lock_keys( %$permissive_peak_calling );
#   
#   my $permissive_peak_file = $permissive_peak_calling->{'peak_file'};
#   my $alignment_name       = $permissive_peak_calling->{'alignment_name'};

  $self->dataflow_output_id(
    {
      execution_plan          => $execution_plan,
      species                 => $species,
      permissive_peak_calling => $permissive_peak_calling,
    }, 
    BRANCH_OUTPUT
  );
}

1;
