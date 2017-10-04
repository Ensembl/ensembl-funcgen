package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRPlanBuilder::SkipIdr;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

sub set_Alignment {
  my $self = shift;
  my $obj  = shift;
  
  return $self->_generic_set('alignment', undef, $obj);
}

sub get_Alignment { 
  return shift->_generic_get('alignment');
}

sub construct {
  my $self = shift;
  my $param = shift;
  my $experiment = $param->{experiment};
  
  $self->set_Alignment([]);
  
  return {
    strategy => Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy->SKIP_IDR,
    type => 'idr',
    name => $experiment->name,
  };
}

1;
