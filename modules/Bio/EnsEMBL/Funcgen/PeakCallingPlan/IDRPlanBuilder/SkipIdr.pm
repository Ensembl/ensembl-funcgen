package Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::SkipIdr;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw( :all );

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

sub set_Alignments {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('alignment', undef, $obj);
}

sub get_Alignments { 
  return shift->_generic_get('alignment');
}

sub set_idr_plan {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('idr_plan', undef, $obj);
}

sub get_idr_plan { 
  return shift->_generic_get('idr_plan');
}

sub construct {
  my $self = shift;
  my $param = shift;
  my $experiment = $param->{experiment};
  
  my $idr_plan = {
    strategy => SKIP_IDR,
    type     => 'idr',
    name     => $experiment->name,
  };
  $self->set_Alignments([]);
  $self->set_idr_plan($idr_plan);
  return; 
}

1;
