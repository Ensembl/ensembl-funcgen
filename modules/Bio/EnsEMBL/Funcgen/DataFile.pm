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

package Bio::EnsEMBL::Funcgen::DataFile;

use strict;

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    data_file_id => 'data_file_id',
    table_id     => 'table_id',
    table_name   => 'table_name',
    path         => 'path',
    file_type    => 'file_type',
    md5sum       => 'md5sum',
  };
}

sub dbID         { return shift->_generic_get_or_set('dbID',         @_); }
sub db           { return shift->_generic_get_or_set('db',           @_); }
sub data_file_id { return shift->_generic_get_or_set('data_file_id', @_); }
sub table_id     { return shift->_generic_get_or_set('table_id',     @_); }
sub table_name   { return shift->_generic_get_or_set('table_name',   @_); }
sub path         { return shift->_generic_get_or_set('path',         @_); }
sub file_type    { return shift->_generic_get_or_set('file_type',    @_); }
sub md5sum       { return shift->_generic_get_or_set('md5sum',       @_); }

1;
