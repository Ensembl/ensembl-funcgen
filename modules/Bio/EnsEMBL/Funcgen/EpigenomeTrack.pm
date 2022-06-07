=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::EpigenomeTrack - Object that represents a peak calling

=head1 SYNOPSIS

  use strict;
  use warnings;
  use Bio::EnsEMBL::Registry;
  use List::Util qw( min );
  use Data::Dumper;

  my $registry = 'Bio::EnsEMBL::Registry';

  $registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org',
      -user => 'anonymous'
  );

  my $epigenome_track_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'EpigenomeTrack');

  my $all_epigenome_tracks = $epigenome_track_adaptor->fetch_all;

  my $number_of_epigenome_tracks_available = @$all_epigenome_tracks;

  print "There are $number_of_epigenome_tracks_available peak callings available for querying:\n";

  # Print the first ten
  my $max_features_to_print = 10;

  for my $i ( 1.. min($max_features_to_print, $number_of_epigenome_tracks_available) ) {

    my $current_epigenome_track = $all_epigenome_tracks->[$i];
    print "  - " . $current_epigenome_track->display_label . "\n";
    
  }

=head1 DESCRIPTION

This object represents a peak calling from a ChIP-seq or other high-throughput
assay. It links to 

  - the set of Peaks that were generated and 
  - the alignment that the peak calling was done on and
  - the peak caller that was used via the analysis.

=cut

package Bio::EnsEMBL::Funcgen::EpigenomeTrack;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
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
   dbID                      => 'dbID',
   adaptor                   => 'adaptor',
   feature_type_id           => 'feature_type_id',
   track_type                => 'track_type',
   epigenome_id              => 'epigenome_id',
   data_file_id              => 'data_file_id'
  };
}

sub dbID            { return shift->_generic_get_or_set('dbID',             @_); }
sub adaptor         { return shift->_generic_get_or_set('adaptor',          @_); }
sub feature_type_id { return shift->_generic_get_or_set('feature_type_id',  @_); }
sub track_type      { return shift->_generic_get_or_set('track_type',       @_); }
sub epigenome_id    { return shift->_generic_get_or_set('epigenome_id',     @_); }
sub data_file_id    { return shift->_generic_get_or_set('data_file_id',     @_); }


=head2 get_FeatureType

  Example    : my $feature_type = $epigenome_track->get_FeatureType;
  Description: Fetches the feature type of the peak calling. This is the
               type of feature the experiment was assaying for.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub get_FeatureType {
  return shift->_generic_fetch('feature_type', 'get_FeatureTypeAdaptor', 'feature_type_id');
}

=head2 get_Epigenome

  Example    : my $epigenome = $epigenome_track->get_Epigenome;
  Description: Gets the epigenome that was used in the assay.
  Returntype : Bio::EnsEMBL::Funcgen::Epigenome
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub get_Epigenome {
    my $self = shift;
    return $self->_generic_fetch('epigenome', 'get_EpigenomeAdaptor', 'epigenome_id');
}

=head2 get_DataFile

  Example     : my $data_file = $epigenome_track->get_DataFile;
  Description : Gets the data file related to the epigenome_track.
  Returntype  : Bio::EnsEMBL::Funcgen::DataFile
  Exception   : None
  Caller      : general
  Status      : Stable

=cut

sub get_DataFile {

  my $self = shift;
  return $self->_generic_fetch('data_file', 'get_DataFileAdaptor', 'data_file_id');

}

=head2 summary_as_hash

  Example       : $summary = $epigenome_track->summary_as_hash;
  Description   : Returns summary in a hash reference.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self   = shift;
  
  my $epigenome         = $self->get_Epigenome;
  my $feature_type      = $self->get_FeatureType;
  my $data_file         = $self->get_DataFile;
  
  
  my $summary = {
    epigenome         => $epigenome     ->summary_as_hash,
    feature_type      => $feature_type  ->summary_as_hash,
    data_file         => $data_file     ->summary_as_hash,
  };
  
  return $summary;
}

1;
