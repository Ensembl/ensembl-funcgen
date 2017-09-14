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

package Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration;

use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    technical_replicate  => 'technical_replicate',
    biological_replicate => 'biological_replicate',
    paired_end_tag       => 'paired_end_tag',
    multiple             => 'multiple',
    read_file_id         => 'read_file_id',
    experiment_id        => 'experiment_id',
  }
}

sub dbID                 { return shift->_generic_get_or_set('dbID',                 @_); }
sub db                   { return shift->_generic_get_or_set('db',                   @_); }
sub biological_replicate { return shift->_generic_get_or_set('biological_replicate', @_); }
sub technical_replicate  { return shift->_generic_get_or_set('technical_replicate',  @_); }
sub paired_end_tag       { return shift->_generic_get_or_set('paired_end_tag',       @_); }
sub multiple             { return shift->_generic_get_or_set('multiple',             @_); }
sub read_file_id         { return shift->_generic_get_or_set('read_file_id',         @_); }
sub experiment_id        { return shift->_generic_get_or_set('experiment_id',        @_); }

sub fetch_ReadFile {

  my $self = shift;
  
  my $read_file_id = $self->read_file_id;
  my $read_file    = $self->db->db->get_ReadFileAdaptor
    ->fetch_by_dbID($read_file_id);

  return $read_file;
}

1;





