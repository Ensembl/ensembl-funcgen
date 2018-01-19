=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration;

use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_get
  _generic_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    technical_replicate  => 'technical_replicate',
    biological_replicate => 'biological_replicate',
    paired_end_tag       => 'paired_end_tag',
    multiple             => 'multiple',
    read_file_id         => '_read_file_id',
    read_file            => 'set_ReadFile',
    experiment_id        => 'experiment_id',
  }
}

sub dbID                 { return shift->_generic_get_or_set('dbID',                 @_); }
sub db                   { return shift->_generic_get_or_set('db',                   @_); }

=head2 biological_replicate

  Example    : my $biological_replicate 
                 = $read_file_experimental_configuration
                   ->biological_replicate;
  Description: Accessor for the biological_replicate number of the ReadFile 
               linked to this object.
  Returntype : Int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub biological_replicate { return shift->_generic_get_or_set('biological_replicate', @_); }

=head2 technical_replicate

  Example    : my $technical_replicate 
                 = $read_file_experimental_configuration
                   ->technical_replicate;
  Description: Accessor for the technical_replicate number of the ReadFile 
               linked to this object.
  Returntype : Int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub technical_replicate  { return shift->_generic_get_or_set('technical_replicate',  @_); }
sub paired_end_tag       { return shift->_generic_get_or_set('paired_end_tag',       @_); }
sub multiple             { return shift->_generic_get_or_set('multiple',             @_); }
sub experiment_id        { return shift->_generic_get_or_set('experiment_id',        @_); }
sub _read_file_id        { return shift->_generic_get_or_set('_read_file_id',        @_); }

=head2 get_ReadFile

  Example    : my $read_file 
                 = $read_file_experimental_configuration
                   ->get_ReadFile;
  Description: Getter for the ReadFile object.
  Returntype : Bio::EnsEMBL::Funcgen::ReadFile
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub get_ReadFile {
  return shift->_generic_get('read_file');
}

=head2 set_ReadFile

  Example    : $read_file_experimental_configuration
                   ->set_ReadFile($read_file);
  Description: Setter for the ReadFile object.
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub set_ReadFile {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('read_file', 'Bio::EnsEMBL::Funcgen::ReadFile', $obj);
}

1;
