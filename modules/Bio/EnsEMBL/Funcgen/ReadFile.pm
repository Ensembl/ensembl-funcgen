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
    dbID           => 'dbID',
    name           => 'name',
    is_paired_end  => 'is_paired_end',
    file_size      => 'file_size',
    read_length    => 'read_length',
    md5sum         => 'md5sum',
    file           => 'file',
    notes          => 'notes',
    analysis_id                          => '_analysis_id',
    analysis                             => 'set_Analysis',
    read_file_experimental_configuration => 'set_ReadFileExperimentalConfiguration',
  };
}

sub dbID          { return shift->_generic_get_or_set('dbID',          @_); }
sub db            { return shift->_generic_get_or_set('db',            @_); }

=head2 name

  Example    : my $name = $peak_calling->name;
  Description: Accessor for the name of the read file.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub name          { return shift->_generic_get_or_set('name',          @_); }
sub is_paired_end { return shift->_generic_get_or_set('is_paired_end', @_); }

=head2 file_size

  Example    : my $file_size = $peak_calling->file_size;
  Description: Accessor for the file_size of the read file.
  Returntype : Int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub file_size     { return shift->_generic_get_or_set('file_size',     @_); }

=head2 read_length

  Example    : my $read_length = $read_file->read_length;
  Description: Accessor for the read_length of the read file.
  Returntype : Int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub read_length   { return shift->_generic_get_or_set('read_length',   @_); }

=head2 md5sum

  Example    : my $md5sum = $read_file->md5sum;
  Description: Accessor for the md5sum of the data file.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub md5sum        { return shift->_generic_get_or_set('md5sum',        @_); }
sub file          { return shift->_generic_get_or_set('file',          @_); }
sub notes         { return shift->_generic_get_or_set('notes',         @_); }
sub _analysis_id  { return shift->_generic_get_or_set('_analysis_id',  @_); }

=head2 get_Analysis

  Example    : my $analysis = $read_file->get_Analysis;
  Description: Getter for the analysis of the read file. This is not an in 
               silico analysis, but the the protocol that was used to 
               generate the reads.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub get_Analysis {
  return shift->_generic_get('analysis');
}

=head2 get_Analysis

  Example    : $read_file->set_Analysis($analysis);
  Description: Setter for the analysis of the read file. The analysis should 
               represent the protocol that was used to generate the reads.
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub set_Analysis {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('analysis', 'Bio::EnsEMBL::Analysis', $obj);
}

=head2 get_ReadFileExperimentalConfiguration

  Example    : my $read_file_experimental_configuration 
                 = $read_file->get_ReadFileExperimentalConfiguration;
  Description: Getter for the ReadFileExperimentalConfiguration object.
  Returntype : Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub get_ReadFileExperimentalConfiguration {
  return shift->_generic_get('read_file_experimental_configuration');
}

=head2 set_ReadFileExperimentalConfiguration

  Example    : $read_file->get_ReadFileExperimentalConfiguration(
                 $read_file_experimental_configuration
               );
  Description: Setter for the ReadFileExperimentalConfiguration object.
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
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
