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

package Bio::EnsEMBL::Funcgen::PhantomPeak;

use strict;

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID         => 'dbID',
    db           => 'db',

    analysis_id         => 'analysis_id',
    alignment_id        => 'alignment_id',
    num_reads           => 'num_reads',
    est_frag_len        => 'est_frag_len',
    est_frag_len_2      => 'est_frag_len_2',
    est_frag_len_3      => 'est_frag_len_3',
    corr_est_frag_len   => 'corr_est_frag_len',
    corr_est_frag_len_2 => 'corr_est_frag_len_2',
    corr_est_frag_len_3 => 'corr_est_frag_len_3',
    phantom_peak        => 'phantom_peak',
    corr_phantom_peak   => 'corr_phantom_peak',
    argmin_corr         => 'argmin_corr',
    min_corr            => 'min_corr',
    nsc                 => 'nsc',
    rsc                 => 'rsc',
    quality_tag         => 'quality_tag'
  };
}

sub dbID         { return shift->_generic_get_or_set('dbID',         @_); }
sub db           { return shift->_generic_get_or_set('db',           @_); }

sub analysis_id         { return shift->_generic_get_or_set(' analysis_id',         @_); }
sub alignment_id        { return shift->_generic_get_or_set(' alignment_id',        @_); }
sub num_reads           { return shift->_generic_get_or_set(' num_reads',           @_); }
sub est_frag_len        { return shift->_generic_get_or_set(' est_frag_len',        @_); }
sub est_frag_len_2      { return shift->_generic_get_or_set(' est_frag_len_2',      @_); }
sub est_frag_len_3      { return shift->_generic_get_or_set(' est_frag_len_3',      @_); }
sub corr_est_frag_len   { return shift->_generic_get_or_set(' corr_est_frag_len',   @_); }
sub corr_est_frag_len_2 { return shift->_generic_get_or_set(' corr_est_frag_len_2', @_); }
sub corr_est_frag_len_3 { return shift->_generic_get_or_set(' corr_est_frag_len_3', @_); }
sub phantom_peak        { return shift->_generic_get_or_set(' phantom_peak',        @_); }
sub corr_phantom_peak   { return shift->_generic_get_or_set(' corr_phantom_peak',   @_); }
sub argmin_corr         { return shift->_generic_get_or_set(' argmin_corr',         @_); }
sub min_corr            { return shift->_generic_get_or_set(' min_corr',            @_); }
sub nsc                 { return shift->_generic_get_or_set(' nsc',                 @_); }
sub rsc                 { return shift->_generic_get_or_set(' rsc',                 @_); }
sub quality_tag         { return shift->_generic_get_or_set(' quality_tag',         @_); }

1;
