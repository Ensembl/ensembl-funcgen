package Bio::EnsEMBL::Funcgen::PeakCallingPlan::SignalFilePlanBuilder;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

# Input

sub set_alignment       { return shift->_generic_set('alignment',       undef, @_); }
sub set_alignment_namer { return shift->_generic_set('alignment_namer', undef, @_); }
sub set_is_control      { return shift->_generic_set('is_control',      undef, @_); }

# Output

sub get_signal_plan { 
  return shift->_generic_get('signal_plan');
}

sub _set_signal_plan {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('signal_plan', undef, $obj);
}

sub _get_alignment       { return shift->_generic_get('alignment')       }
sub _get_alignment_namer { return shift->_generic_get('alignment_namer') }
sub _get_is_control      { return shift->_generic_get('is_control')      }

sub construct {
  my $self  = shift;
  
  my $alignment_namer = $self->_get_alignment_namer;
  
  my $is_control = $self->_get_is_control;
  
  my $is_valid = ($is_control eq TRUE) || ($is_control eq FALSE);
  
  if (! $is_valid) {
    use Carp;
    confess("Invalid value for is_control! $is_control!");
  }
  
   my $bigwig = {
    input    => $self->_get_alignment,
    type     => SIGNAL_EXPERIMENT,
    name     => $alignment_namer->base_name_no_duplicates,
    analysis => CONVERT_BAM_TO_BIGWIG_ANALYSIS,
    is_control  => $is_control,
    output   => {
      real   => $alignment_namer->bigwig_file_no_duplicates,
      stored => $alignment_namer->bigwig_file_no_duplicates_stored,
      format => BIGWIG_FORMAT
    }
  };
  $self->_set_signal_plan($bigwig);
  return;
}

1;
