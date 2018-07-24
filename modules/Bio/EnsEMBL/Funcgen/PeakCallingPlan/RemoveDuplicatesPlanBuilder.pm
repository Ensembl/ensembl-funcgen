package Bio::EnsEMBL::Funcgen::PeakCallingPlan::RemoveDuplicatesPlanBuilder;

use strict;
use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

# Input

sub set_input         { return shift->_generic_set('input',         undef, @_); }
sub set_name          { return shift->_generic_set('name',          undef, @_); }
sub set_output_real   { return shift->_generic_set('output_real',   undef, @_); }
sub set_output_stored { return shift->_generic_set('output_stored', undef, @_); }
sub set_output_format { return shift->_generic_set('output_format', undef, @_); }
sub set_is_control    { return shift->_generic_set('is_control',    undef, @_); }
sub set_experiment    { return shift->_generic_set('experiment',    undef, @_); }

# Output

sub get_plan { return shift->_generic_get('plan'); }

sub _get_input         { return shift->_generic_get('input'); }
sub _get_name          { return shift->_generic_get('name'); }
sub _get_output_real   { return shift->_generic_get('output_real'); }
sub _get_output_stored { return shift->_generic_get('output_stored'); }
sub _get_output_format { return shift->_generic_get('output_format'); }
sub _get_is_control    { return shift->_generic_get('is_control'); }
sub _get_experiment    { return shift->_generic_get('experiment'); }

sub _set_plan { return shift->_generic_set('plan', undef, @_); }

sub construct {
  my $self = shift;
  
  my $is_control = $self->_get_is_control;
  
  my $is_valid = ($is_control eq TRUE) || ($is_control eq FALSE);
  
  if (! $is_valid) {
    use Carp;
    confess("Invalid value for is_control! $is_control!");
  }

  my $remove_duplicates_plan = {
    input           => $self->_get_input,
    name            => $self->_get_name,
    type            => ALIGNMENT_TYPE,
    analysis        => REMOVE_DUPLICATES_ANALYSIS,
    is_control      => $is_control,
    from_experiment => $self->_get_experiment,
    output => {
      real   => $self->_get_output_real,
      stored => $self->_get_output_stored,
      format => $self->_get_output_format,
    },
  };
  $self->_set_plan($remove_duplicates_plan);
  return;
}

1;
