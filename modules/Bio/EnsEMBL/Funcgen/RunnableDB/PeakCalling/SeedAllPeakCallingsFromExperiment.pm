package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllPeakCallingsFromExperiment;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_OUTPUT => 2,
};

sub run {

  my $self           = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');
  
  my $peak_calling_name = $execution_plan->{call_peaks}->{name};
  
  $self->dataflow_output_id( 
    {
      'peak_calling'  => $peak_calling_name,
      'species'       => $species,
    }, 
    BRANCH_OUTPUT
  );
}

1;
