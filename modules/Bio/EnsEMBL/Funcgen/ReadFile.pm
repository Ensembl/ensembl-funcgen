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
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
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
    dbID            => 'dbID',
    adaptor         => 'adaptor',
    name            => 'name',
    is_paired_end   => 'is_paired_end',
    file_size       => 'file_size',
    number_of_reads => 'number_of_reads',
    read_length     => 'read_length',
    md5sum          => 'md5sum',
    file            => 'file',
    notes           => 'notes',
    analysis_id     => '_analysis_id',
    analysis        => 'set_Analysis',
  };
}

sub dbID          { return shift->_generic_get_or_set('dbID',          @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}

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
sub md5sum          { return shift->_generic_get_or_set('md5sum',          @_); }
sub file            { return shift->_generic_get_or_set('file',            @_); }
sub number_of_reads { return shift->_generic_get_or_set('number_of_reads', @_); }
sub notes           { return shift->_generic_get_or_set('notes',           @_); }
sub _analysis_id    { return shift->_generic_get_or_set('_analysis_id',    @_); }

sub paired_end_tag {
  my $self = shift;
  if (! $self->is_paired_end) {
    throw("Not a paired end read file!");
  }
  my $read_file_experimental_configuration = $self->get_ReadFileExperimentalConfiguration;
  my $paired_end_tag = $read_file_experimental_configuration->paired_end_tag;
  return $paired_end_tag;
}

sub get_mate_ReadFile {
  my $self = shift;
  
  if (! $self->is_paired_end) {
    throw("Not a paired end read file!");
  }
  
  my $paired_end_tag      =     $self->paired_end_tag;
  my $mate_paired_end_tag = 3 - $paired_end_tag;
  
  my $read_file_experimental_configuration = $self->get_ReadFileExperimentalConfiguration;

  my $pairs_read_file_experimental_configuration 
    = Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration->new(
      -technical_replicate  => $read_file_experimental_configuration->technical_replicate,
      -biological_replicate => $read_file_experimental_configuration->biological_replicate,
      -paired_end_tag       => $mate_paired_end_tag,
      -multiple             => $read_file_experimental_configuration->multiple,
      -experiment_id        => $read_file_experimental_configuration->experiment_id,
    );
  
  my $read_file_adaptor 
    = $self->adaptor->db->get_ReadFileAdaptor;
  
  my $read_file_mate
    = $read_file_adaptor
      ->fetch_by_ReadFileExperimentalConfiguration(
        $pairs_read_file_experimental_configuration
      );
  
  if (! defined $read_file_mate) {
    throw("Couldn't find pair!");
  }
  return $read_file_mate;
}

sub get_FastQC {

  my $self         = shift;
  
  my $fastqc_adaptor = $self->adaptor->db->get_FastQCAdaptor;
  if (! defined $fastqc_adaptor) {
    throw("Couldn't get an FastQCAdaptor!");
  }
  my $fastqc = $fastqc_adaptor->fetch_by_ReadFile($self);
  return $fastqc;
}

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

  my $self = shift;

  my $read_file_experimental_configuration_adaptor
    = $self->adaptor->db->get_ReadFileExperimentalConfigurationAdaptor;

  my $read_file_experimental_configurations 
    = $read_file_experimental_configuration_adaptor
      ->fetch_all_by_read_file_id(
        $self->dbID
      );
  return $read_file_experimental_configurations->[0];
}

=head2 set_ReadFileExperimentalConfiguration

  Example    : $read_file->fetch_ReadFileExperimentalConfiguration(
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

=head2 summary_as_hash

  Example       : $summary = $read_file->summary_as_hash;
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
  $suppress_link = '' if (! defined $suppress_link);
  
  my $summary = {
    name => $self->name,
  };
  
  if ($suppress_link ne 'fastqc') {
    my $fastqc = $self->get_FastQC;
    if (defined $fastqc) {
      $summary->{'fastqc'} = $fastqc->summary_as_hash('read_file');
    }
  }
  
  return $summary;
}


1;
