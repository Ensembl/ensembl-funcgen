=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Funcgen::Alignment;

use strict;
use warnings;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );

sub _constructor_parameters {
  return {
    dbID           => 'dbID',
    adaptor        => 'adaptor',
    name           => 'name',
    analysis_id    => 'analysis_id',
    bam_file_id    => 'bam_file_id',
    bigwig_file_id => 'bigwig_file_id',
    read_file_ids  => 'read_file_ids',
    experiment_id  => 'experiment_id',
    has_duplicates => 'has_duplicates',
    is_control     => 'is_control',
    to_gender      => 'to_gender',
    is_complete    => 'is_complete',
    
    source_alignment_id       => 'source_alignment_id',
    deduplicated_alignment_id => 'deduplicated_alignment_id',
    
  };
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_set
  _generic_get
  _generic_fetch
);

sub dbID            { return shift->_generic_get_or_set('dbID',           @_); }
sub adaptor         { return shift->_generic_get_or_set('adaptor',        @_); }

=head2 name

  Example    : my $name = $alignment->name;
  Description: Accessor for the name of the alignment.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub name           { return shift->_generic_get_or_set('name',                    @_); }
sub analysis_id    { return shift->_generic_get_or_set('analysis_id',             @_); }
sub bam_file_id    { return shift->_generic_get_or_set('bam_file_id',             @_); }
sub bigwig_file_id { return shift->_generic_get_or_set('bigwig_file_id',          @_); }
sub experiment_id  { return shift->_generic_get_or_set('experiment_id',           @_); }
sub read_file_ids  { return shift->_generic_get_or_set('read_file_ids',           @_); }
sub has_duplicates { return shift->_generic_get_or_set('has_duplicates',          @_); }
sub is_control     { return shift->_generic_get_or_set('is_control',              @_); }
sub to_gender      { return shift->_generic_get_or_set('to_gender',               @_); }
sub is_complete    { return shift->_generic_get_or_set('is_complete',             @_); }

sub source_alignment_id       { return shift->_generic_get_or_set('source_alignment_id',       @_); }
sub deduplicated_alignment_id { return shift->_generic_get_or_set('deduplicated_alignment_id', @_); }

=head2 get_Analysis

  Example    : my $analysis = $alignment->get_Analysis;
  Description: Gets the analysis of the alignment. This is the analysis
               representing the aligner that was used to align the 
               reads.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub get_Analysis {
  return shift->_generic_fetch('analysis', 'get_AnalysisAdaptor', 'analysis_id');
}

sub set_Analysis {
  my $self = shift;
  my $obj  = shift;
  return shift->_generic_set('analysis', 'Bio::EnsEMBL::Analysis', $obj);
}

sub _fetch_DataFile {

  my $self         = shift;
  my $data_file_id = shift;
  
  my $data_file_adaptor = $self->adaptor->db->get_DataFileAdaptor;
  if (! defined $data_file_adaptor) {
    throw("Couldn't get a DataFileAdaptor!");
  }
  my $data_file = $data_file_adaptor->fetch_by_dbID($data_file_id);
  return $data_file;
}

sub get_all_ReadFileExperimentalConfigurations {

  my $self = shift;
  my @read_file_experimental_configurations;
  
  my $read_files  = $self->get_all_ReadFiles;
  foreach my $read_file (@$read_files) {
  
    my $read_file_experimental_configuration 
      = $read_file->get_ReadFileExperimentalConfiguration;
    
    push @read_file_experimental_configurations,
      $read_file_experimental_configuration;
  }
  return \@read_file_experimental_configurations;
}

sub _delete_bam_file_from_db {

  my $self = shift;
  
  $self->_generic_get_or_set('bam_file_id', undef, 1);
  $self->adaptor->update($self);
}

sub get_all_deduplicated_replicate_Alignments {
  my $self = shift;
  
  my $alignment_adaptor = $self->adaptor;
  if (! defined $alignment_adaptor) {
    throw("Couldn't get a AlignmentAdaptor!");
  }
  return $alignment_adaptor->fetch_all_deduplicated_replicates_by_Alignment($self);
}

sub get_source_Alignment {

  my $self         = shift;
  
  my $alignment_adaptor = $self->adaptor;
  if (! defined $alignment_adaptor) {
    throw("Couldn't get a AlignmentAdaptor!");
  }
  my $alignment = $alignment_adaptor->fetch_by_dbID($self->source_alignment_id);
  return $alignment;
}

sub get_Chance_by_control_Alignment {
  my $self = shift;
  my $control_alignment = shift;
  
  my $chance_adaptor = $self->adaptor->db->get_ChanceAdaptor;
  if (! defined $chance_adaptor) {
    throw("Couldn't get an ChanceAdaptor!");
  }
  my $chance = $chance_adaptor
    ->fetch_by_signal_control_Alignments(
      $self, 
      $control_alignment
    );
  return $chance;
}

sub get_Experiment {

  my $self = shift;
  
  my $experiment_adaptor = $self->adaptor->db->get_ExperimentAdaptor;
  if (! defined $experiment_adaptor) {
    throw("Couldn't get a ExperimentAdaptor!");
  }
  my $experiment = $experiment_adaptor->fetch_by_dbID($self->experiment_id);
  return $experiment;
}

sub get_PhantomPeak {
  my $self = shift;

  my $phantom_peak_adaptor = $self->adaptor->db->get_PhantomPeakAdaptor;
  if (! defined $phantom_peak_adaptor) {
    throw("Couldn't get an PhantomPeakAdaptor!");
  }
  my $phantom_peak = $phantom_peak_adaptor
    ->fetch_by_Alignment($self);
  return $phantom_peak;
}

=head2 has_bam_DataFile

  Example    : my $bam_DataFile = undef;
               if ($alignment->has_bam_DataFile) {
                 $bam_DataFile = $alignment->get_bam_DataFile;
               } else {
                 warn "No bam file available!";
               }
  Description: Indicates whether a bam file is available for this alignment.
  Returntype : Boolean
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub has_bam_DataFile {

  my $self = shift;
  return defined $self->bam_file_id;
}

=head2 has_bigwig_DataFile

  Example    : my $bigwig_DataFile = undef;
               if ($alignment->has_bigwig_DataFile) {
                 $bigwig_DataFile = $alignment->get_bigwig_DataFile;
               } else {
                 warn "No bigwig file available!";
               }
  Description: Indicates whether a bam file is available for this alignment.
  Returntype : Boolean
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub has_bigwig_DataFile {

  my $self = shift;
  return defined $self->bigwig_file_id;
}

=head2 get_bam_DataFile

  Example    : my $bam_DataFile = undef;
               if ($alignment->has_bam_DataFile) {
                 $bam_DataFile = $alignment->get_bam_DataFile;
               } else {
                 warn "No bam file available!";
               }
  Description: Gets the data file object representing the bam file for
               this alignment.
               reads.
  Returntype : Bio::EnsEMBL::Funcgen::DataFile
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub get_bam_DataFile {

  my $self        = shift;
  my $bam_file_id = $self->bam_file_id;
  
  return $self->_fetch_DataFile($bam_file_id);
}

=head2 get_bigwig_DataFile

  Example    : my $bigwig_DataFile = undef;
               if ($alignment->has_bigwig_DataFile) {
                 $bigwig_DataFile = $alignment->get_bigwig_DataFile;
               } else {
                 warn "No bigwig file available!";
               }
  Description: Gets the data file object representing the signal file in
               bigwig format for this alignment.
               reads.
  Returntype : Bio::EnsEMBL::Funcgen::DataFile
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub get_bigwig_DataFile {

  my $self           = shift;
  my $bigwig_file_id = $self->bigwig_file_id;
  
  return $self->_fetch_DataFile($bigwig_file_id);
}

=head2 get_all_ReadFiles

  Example    : my $read_files = $alignment->get_all_ReadFiles;
               print @$read_files ." read files were used to generate this alignment.";
               foreach my $current_readfile (@$read_files) {
                 print "  - " . $current_readfile->name . "\n";
               }
  Description: Gets all read file objects representing the read files that
               were used to generate this alignment.
  Returntype : ArrayRef[Bio::EnsEMBL::Funcgen::ReadFile]
  Caller     : general
  Status     : Stable

=cut

sub get_all_ReadFiles {

  my $self = shift;
  
  my $read_file_id = $self->read_file_ids;
  my $read_file_adaptor = $self->adaptor->db->get_ReadFileAdaptor;
  
  if (! defined $read_file_adaptor) {
    throw("Couldn't get a ReadFileAdaptor!");
  }
  
  my @all_read_files;
  
  foreach my $current_read_file_id (@$read_file_id) {
    my $current_read_file = $read_file_adaptor->fetch_by_dbID($current_read_file_id);
    push @all_read_files, $current_read_file;
  }
  return \@all_read_files;
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
  $suppress_link = '' if (! defined $suppress_link);
  
  my $bam_file    = $self->get_bam_DataFile;
  my $bigwig_file = $self->get_bigwig_DataFile;
  
  my $summary = {
    name           => $self->name,
    has_duplicates => $self->has_duplicates,
    is_control     => $self->is_control,
    to_gender      => $self->to_gender,
    is_complete    => $self->is_complete,
  };
  
  if ($bam_file) {
    $summary->{'bam_file'} = $bam_file->summary_as_hash('alignment');
  }
  if ($bigwig_file) {
    $summary->{'bigwig_file'} = $bigwig_file->summary_as_hash('alignment');
  }
  if ($suppress_link ne 'phantom_peak') {
    my $phantom_peak = $self->get_PhantomPeak;
    if (defined $phantom_peak) {
      $summary->{'phantom_peak'} = $phantom_peak->summary_as_hash('alignment');
    }
  }
  if ($suppress_link ne 'read_file') {
    my $read_files = $self->get_all_ReadFiles;
    $summary->{'read_files'} = [ map { $_->summary_as_hash } @$read_files ];
  }
  
  return $summary;
}

1;

