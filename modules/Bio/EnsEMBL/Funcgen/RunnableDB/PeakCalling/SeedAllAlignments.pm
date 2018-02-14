package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllAlignments;

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
  
  my $alignment_to_plan_hash = $execution_plan->{alignment};
  my @alignment_names = keys %$alignment_to_plan_hash;
  
  foreach my $alignment_name (@alignment_names) {
    $self->dataflow_output_id( 
      {
        'alignment' => $alignment_name,
        'species'   => $species,
      }, 
      BRANCH_OUTPUT
    );
  }
}

1;
