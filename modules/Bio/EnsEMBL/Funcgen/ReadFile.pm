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

package Bio::EnsEMBL::Funcgen::ReadFile;

use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_set
  _generic_get
  _generic_fetch
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    name           => 'name',
    is_paired_end  => 'is_paired_end',
    paired_with    => 'paired_with',
    file_size      => 'file_size',
    read_length    => 'read_length',
    md5sum         => 'md5sum',
    file           => 'file',
    notes          => 'notes',
    analysis                             => 'set_Analysis',
    read_file_experimental_configuration => 'set_ReadFileExperimentalConfiguration',
  };
}

sub dbID          { return shift->_generic_get_or_set('dbID',          @_); }
sub db            { return shift->_generic_get_or_set('db',            @_); }
sub name          { return shift->_generic_get_or_set('name',          @_); }
sub is_paired_end { return shift->_generic_get_or_set('is_paired_end', @_); }
sub paired_with   { return shift->_generic_get_or_set('paired_with',   @_); }
sub file_size     { return shift->_generic_get_or_set('file_size',     @_); }
sub read_length   { return shift->_generic_get_or_set('read_length',   @_); }
sub md5sum        { return shift->_generic_get_or_set('md5sum',        @_); }
sub file          { return shift->_generic_get_or_set('file',          @_); }
sub notes         { return shift->_generic_get_or_set('notes',         @_); }

sub _get_methods {
  return [
    {
      method_name => 'get_Analysis',
      hash_key    => 'analysis',
    },
    {
      method_name => 'get_ReadFileExperimentalConfiguration',
      hash_key    => 'read_file_experimental_configuration',
    },
  ]
}

sub get_Analysis {
  return shift->_generic_get('analysis');
}

sub set_Analysis {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('analysis', 'Bio::EnsEMBL::Analysis', $obj);
}

sub get_ReadFileExperimentalConfiguration {
  return shift->_generic_get('read_file_experimental_configuration');
}

sub set_ReadFileExperimentalConfiguration {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set(
    'read_file_experimental_configuration', 
    'Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration', 
    $obj
  );
}


1;
