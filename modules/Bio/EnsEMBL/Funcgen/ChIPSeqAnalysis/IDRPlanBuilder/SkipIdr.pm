package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRPlanBuilder::SkipIdr;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub construct {
  my $self = shift;
  return {
    type => Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy->SKIP_IDR
  };
}

1;
