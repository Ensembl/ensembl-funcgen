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

package Bio::EnsEMBL::Funcgen::FastQC;

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
    adaptor      => 'adaptor',

    read_file_id                     => 'read_file_id',
    basic_statistics                 => 'basic_statistics',
    per_base_sequence_quality        => 'per_base_sequence_quality',
    per_tile_sequence_quality        => 'per_tile_sequence_quality',
    per_sequence_quality_scores      => 'per_sequence_quality_scores',
    per_base_sequence_content        => 'per_base_sequence_content',
    per_sequence_gc_content          => 'per_sequence_gc_content',
    per_base_n_content               => 'per_base_n_content',
    sequence_length_distribution     => 'sequence_length_distribution',
    sequence_duplication_levels      => 'sequence_duplication_levels',
    overrepresented_sequences        => 'overrepresented_sequences',
    adapter_content                  => 'adapter_content',
    kmer_content                     => 'kmer_content',
    run_failed                       => 'run_failed',
    error_message                    => 'error_message',
  };
}

sub dbID         { return shift->_generic_get_or_set('dbID',         @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub read_file_id                 { return shift->_generic_get_or_set('read_file_id',                 @_); }
sub basic_statistics             { return shift->_generic_get_or_set('basic_statistics',             @_); }
sub per_base_sequence_quality    { return shift->_generic_get_or_set('per_base_sequence_quality',    @_); }
sub per_tile_sequence_quality    { return shift->_generic_get_or_set('per_tile_sequence_quality',    @_); }
sub per_sequence_quality_scores  { return shift->_generic_get_or_set('per_sequence_quality_scores',  @_); }
sub per_base_sequence_content    { return shift->_generic_get_or_set('per_base_sequence_content',    @_); }
sub per_sequence_gc_content      { return shift->_generic_get_or_set('per_sequence_gc_content',      @_); }
sub per_base_n_content           { return shift->_generic_get_or_set('per_base_n_content',           @_); }
sub sequence_length_distribution { return shift->_generic_get_or_set('sequence_length_distribution', @_); }
sub sequence_duplication_levels  { return shift->_generic_get_or_set('sequence_duplication_levels',  @_); }
sub overrepresented_sequences    { return shift->_generic_get_or_set('overrepresented_sequences',    @_); }
sub adapter_content              { return shift->_generic_get_or_set('adapter_content',              @_); }
sub kmer_content                 { return shift->_generic_get_or_set('kmer_content',                 @_); }
sub run_failed                   { return shift->_generic_get_or_set('run_failed',                   @_); }
sub error_message                { return shift->_generic_get_or_set('error_message',                @_); }

=head2 summary_as_hash

  Example       : $summary = $fastqc->summary_as_hash;
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
    basic_statistics                 => $self->basic_statistics,
    per_base_sequence_quality        => $self->per_base_sequence_quality,
    per_tile_sequence_quality        => $self->per_tile_sequence_quality,
    per_sequence_quality_scores      => $self->per_sequence_quality_scores,
    per_base_sequence_content        => $self->per_base_sequence_content,
    per_sequence_gc_content          => $self->per_sequence_gc_content,
    per_base_n_content               => $self->per_base_n_content,
    sequence_length_distribution     => $self->sequence_length_distribution,
    sequence_duplication_levels      => $self->sequence_duplication_levels,
    overrepresented_sequences        => $self->overrepresented_sequences,
    adapter_content                  => $self->adapter_content,
    kmer_content                     => $self->kmer_content
  };

#   if ($suppress_link ne 'read_file') {
#     my $fastqc = $self->fetch_FastQC;
#     $summary->{'fastqc'} = $fastqc->summary_as_hash('read_file');
#   }

  return $summary;
}

1;
