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
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID                => 'dbID',
    adaptor             => 'adaptor',
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
    quality_tag         => 'quality_tag',
    run_failed          => 'run_failed',
    error_message       => 'error_message',
  };
}

sub dbID         { return shift->_generic_get_or_set('dbID',         @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub analysis_id         { return shift->_generic_get_or_set('analysis_id',         @_); }
sub alignment_id        { return shift->_generic_get_or_set('alignment_id',        @_); }
sub num_reads           { return shift->_generic_get_or_set('num_reads',           @_); }
sub est_frag_len        { return shift->_generic_get_or_set('est_frag_len',        @_); }
sub est_frag_len_2      { return shift->_generic_get_or_set('est_frag_len_2',      @_); }
sub est_frag_len_3      { return shift->_generic_get_or_set('est_frag_len_3',      @_); }
sub corr_est_frag_len   { return shift->_generic_get_or_set('corr_est_frag_len',   @_); }
sub corr_est_frag_len_2 { return shift->_generic_get_or_set('corr_est_frag_len_2', @_); }
sub corr_est_frag_len_3 { return shift->_generic_get_or_set('corr_est_frag_len_3', @_); }
sub phantom_peak        { return shift->_generic_get_or_set('phantom_peak',        @_); }
sub corr_phantom_peak   { return shift->_generic_get_or_set('corr_phantom_peak',   @_); }
sub argmin_corr         { return shift->_generic_get_or_set('argmin_corr',         @_); }
sub min_corr            { return shift->_generic_get_or_set('min_corr',            @_); }
sub nsc                 { return shift->_generic_get_or_set('nsc',                 @_); }
sub rsc                 { return shift->_generic_get_or_set('rsc',                 @_); }
sub quality_tag         { return shift->_generic_get_or_set('quality_tag',         @_); }

sub run_failed          { return shift->_generic_get_or_set('run_failed',          @_); }
sub error_message       { return shift->_generic_get_or_set('error_message',       @_); }

sub get_Alignment {

  my $self = shift;
  
  my $alignment_adaptor = $self->adaptor->db->get_AlignmentAdaptor;
  if (! defined $alignment_adaptor) {
    throw("Couldn't get an AlignmentAdaptor!");
  }
  my $alignment = $alignment_adaptor->fetch_by_dbID($self->alignment_id);
  return $alignment;
}

=head2 summary_as_hash

  Example       : $summary = $peak_calling->summary_as_hash;
  Description   : Returns summary in a hash reference.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self   = shift;
  
  # Optional parameter to avoid infinite recursions when two objects 
  # reference each other.
  #
  my $suppress_link = shift;

  
  my $summary = {
    num_reads    => $self->num_reads,
    quality_tag  => $self->quality_tag,
    rsc          => $self->rsc,
    nsc          => $self->nsc,
    phantom_peak => $self->phantom_peak,
    est_frag_len => $self->est_frag_len,
  };
  
  if ($suppress_link ne 'alignment') {
    my $alignment  = $self->get_Alignment;
    $summary->{'alignment'} = $alignment->summary_as_hash('phantom_peak');
  }
  
  return $summary;
}

1;
