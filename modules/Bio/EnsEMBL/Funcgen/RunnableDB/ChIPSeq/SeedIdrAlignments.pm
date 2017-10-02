package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SeedIdrAlignments;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_ALIGN => 2,
};

sub run {

  my $self = shift;
  
  my $species = $self->param_required('species');
  my $run_idr = $self->param_required('run_idr');
  
  my $alignment_replicates = $run_idr->{alignment_replicates};
  
  for my $alignment_replicate (@$alignment_replicates) {
  
    my $job = {
      plan    => $alignment_replicate,
      species => $species,
    };
    $self->dataflow_output_id( 
      $job,
      BRANCH_ALIGN
    );
  }
}

1;
