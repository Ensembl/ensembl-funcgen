package Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentPlanFactory;

use strict;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    'names_of_reads_to_merge' => 'names_of_reads_to_merge',
    'description'             => 'description',
    'name'                    => 'name',
    'to_gender'               => 'to_gender',
    'to_assembly'             => 'to_assembly',
    'ensembl_analysis'        => 'ensembl_analysis',
    'output_real'             => 'output_real',
    'output_stored'           => 'output_stored',
    'output_format'           => 'output_format',
    'is_control'              => 'is_control',
    'from_experiment'         => 'from_experiment',
    'has_all_reads'           => 'has_all_reads',
  };
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub names_of_reads_to_merge { return shift->_generic_get_or_set('names_of_reads_to_merge', @_); }
sub description             { return shift->_generic_get_or_set('description',             @_); }
sub name                    { return shift->_generic_get_or_set('name',                    @_); }
sub to_gender               { return shift->_generic_get_or_set('to_gender',               @_); }
sub to_assembly             { return shift->_generic_get_or_set('to_assembly',             @_); }
sub ensembl_analysis        { return shift->_generic_get_or_set('ensembl_analysis',        @_); }
sub output_real             { return shift->_generic_get_or_set('output_real',             @_); }
sub output_stored           { return shift->_generic_get_or_set('output_stored',           @_); }
sub output_format           { return shift->_generic_get_or_set('output_format',           @_); }
sub is_control              { return shift->_generic_get_or_set('is_control',              @_); }
sub from_experiment         { return shift->_generic_get_or_set('from_experiment',         @_); }
sub has_all_reads           { return shift->_generic_get_or_set('has_all_reads',           @_); }

sub product {
  my $self = shift;

  my $is_control = $self->is_control;
  
  my $align_plan = {
    input => {
      read_files => $self->names_of_reads_to_merge,
    },
    description => $self->description,
    name        => $self->name,
    to_gender   => $self->to_gender,
    to_assembly => $self->to_assembly,
    type        => ALIGNMENT_TYPE,
    analysis    => ALIGNMENT_ANALYSIS,
    ensembl_analysis => $self->ensembl_analysis,
    from_experiment => $self->from_experiment,
    is_control    => $is_control,
    is_complete => $self->has_all_reads,
    output => {
      real   => $self->output_real,
      stored => $self->output_stored,
      format => $self->output_format,
    },
  };
  return $align_plan;
}

1;
