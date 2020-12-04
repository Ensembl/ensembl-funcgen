=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Chance;

use strict;
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID         => 'dbID',
    adaptor           => 'adaptor',
    chance_id                                                       => 'chance_id',
    signal_alignment_id                                             => 'signal_alignment_id',
    control_alignment_id                                            => 'control_alignment_id',
    analysis_id                                                     => 'analysis_id',
    p                                                               => 'p',
    'q'                                                             => 'q',
    divergence                                                      => 'divergence',
    z_score                                                         => 'z_score',
    percent_genome_enriched                                         => 'percent_genome_enriched',
    input_scaling_factor                                            => 'input_scaling_factor',
    differential_percentage_enrichment                              => 'differential_percentage_enrichment',
    control_enrichment_stronger_than_chip_at_bin                    => 'control_enrichment_stronger_than_chip_at_bin',
    first_nonzero_bin_at                                            => 'first_nonzero_bin_at',
    pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome => 'pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome',
    path                                                            => 'path',
    run_failed          => 'run_failed',
    error_message       => 'error_message',
  };
}

sub dbID {return shift->_generic_get_or_set('dbID', @_);}
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub chance_id                                                       { return shift->_generic_get_or_set('chance_id',               @_); }
sub signal_alignment_id                                             { return shift->_generic_get_or_set('signal_alignment_id',     @_); }
sub control_alignment_id                                            { return shift->_generic_get_or_set('control_alignment_id',    @_); }
sub analysis_id                                                     { return shift->_generic_get_or_set('analysis_id',             @_); }
sub p                                                               { return shift->_generic_get_or_set('p',                       @_); }
sub q                                                               { return shift->_generic_get_or_set('q',                       @_); }
sub divergence                                                      { return shift->_generic_get_or_set('divergence',              @_); }
sub z_score                                                         { return shift->_generic_get_or_set('z_score',                 @_); }
sub percent_genome_enriched                                         { return shift->_generic_get_or_set('percent_genome_enriched', @_); }
sub input_scaling_factor                                            { return shift->_generic_get_or_set('input_scaling_factor',    @_); }
sub differential_percentage_enrichment                              { return shift->_generic_get_or_set('differential_percentage_enrichment',   @_); }
sub control_enrichment_stronger_than_chip_at_bin                    { return shift->_generic_get_or_set('control_enrichment_stronger_than_chip_at_bin',   @_); }
sub first_nonzero_bin_at                                            { return shift->_generic_get_or_set('first_nonzero_bin_at',   @_); }
sub pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome { return shift->_generic_get_or_set('pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome',   @_); }
sub path                                                            { return shift->_generic_get_or_set('path',   @_); }

sub run_failed          { return shift->_generic_get_or_set('run_failed',          @_); }
sub error_message       { return shift->_generic_get_or_set('error_message',       @_); }

=head2 summary_as_hash

  Example       : $summary = $chance->summary_as_hash;
  Description   : Returns summary in a hash reference.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self   = shift;
  
  return {
    'p' => $self->p,
    'q' => $self->q,
    'divergence' => $self->divergence,
    'z_score' => $self->z_score,
    'percent_genome_enriched' => $self->percent_genome_enriched,
    'input_scaling_factor' => $self->input_scaling_factor,
    'differential_percentage_enrichment' => $self->differential_percentage_enrichment,
    'control_enrichment_stronger_than_chip_at_bin' => $self->control_enrichment_stronger_than_chip_at_bin,
  };
}


1;
