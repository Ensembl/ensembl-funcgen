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

Bio::EnsEMBL::Funcgen::RegulatoryBuild - A module to represent a Regulatory Build.

=head1 SYNOPSIS

  use strict;
  use warnings;
  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
  );

  my $regulatory_build_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryBuild');
  my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
  my $epigenomes_in_regulatory_build = $regulatory_build->get_all_Epigenomes;

  print "There are " . @$epigenomes_in_regulatory_build . " epigenomes in the regulatory build:\n";
  print join "\n", map { '  - ' . $_->display_label } @$epigenomes_in_regulatory_build;
  print  "\n";

  print "Name: "                                   . $regulatory_build->name                                      . "\n";
  print "Version: "                                . $regulatory_build->version                                   . "\n";
  print "Initial release date: "                   . $regulatory_build->initial_release_date                      . "\n";
  print "Last annotation update: "                 . $regulatory_build->last_annotation_update                    . "\n";
  print "feature type id: "                        . $regulatory_build->feature_type_id                           . "\n";
  print "Analysis id: "                            . $regulatory_build->analysis_id                               . "\n";
  print "Is current: "                             . $regulatory_build->is_current                                . "\n";
  print "Stable id of sample regulatory feature: " . $regulatory_build->fetch_sample_RegulatoryFeature->stable_id . "\n";

=head1 DESCRIPTION

  This is a class represents a regulatory build.

=cut

package Bio::EnsEMBL::Funcgen::RegulatoryBuild;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use base 'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality';

=head2 new
  Arg [-db] :
  Arg [-dbID] :
  Arg [-name] :
        name of regulatory build
  Arg [-version] :
        version of regulatory build
  Arg [-initial_release_date]   :
        date of the initial release
  Arg [-last_annotation_update] :
        date of the last update
  Arg [-feature_type_id] :
        the feature type id
  Arg [-analysis_id] :
        the analysis id
  Arg [-is_current] :
        is current release

  Example      : my $regulatory_build = Bio::EnsEMBL::Funcgen::RegulatoryBuild->new(
                            -name,                  => 'The Ensembl Regulatory Build', 
                            -version                => '14',
                            -initial_release_date   => '2016-6',
                            -last_annotation_update => '2016-6',
                            -feature_type_id        => 19,
                            -analysis_id            => 16,
                            -is_current             => 1
                          );
  Description : Constructor method for Regulatory Build class
  Returntype  : Bio::EnsEMBL::Funcgen::RegulatoryBuild
  Exceptions  : None 
  Caller      : General
  Status      : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = bless {}, $class;

  my @field = qw(
    db
    dbID
    name
    version
    initial_release_date
    last_annotation_update
    feature_type_id
    feature_type
    analysis_id
    analysis
    is_current
    sample_regulatory_feature_id
  );
  
  my (
    $db,
    $dbID,
    $name,
    $version,
    $initial_release_date,
    $last_annotation_update,
    $feature_type_id,
    $feature_type,
    $analysis_id,
    $analysis,
    $is_current,
    $sample_regulatory_feature_id,
  )
    = rearrange([ @field ], @_);

  #
  # Analysis and FeatureType can be passed in either as objects or by their 
  # dbIDs. The RegulatoryFeatureAdaptor uses database ids only, the objects
  # are lazy loaded as needed.
  #
  # The regulatory build script uses proper objects to create the regulatory
  # build object.
  #
  if (defined $analysis) {
    $self->set_Analysis($analysis);
  }
  if (defined $feature_type) {
    $self->set_FeatureType($feature_type);
  }
  
  $self->feature_type_id($feature_type_id);
  $self->analysis_id($analysis_id);

  $self->db($db);
  $self->dbID($dbID);
  $self->name($name);
  $self->version($version);
  $self->initial_release_date($initial_release_date);
  $self->last_annotation_update($last_annotation_update);
  $self->is_current($is_current);
  $self->sample_regulatory_feature_id($sample_regulatory_feature_id);
  

  return $self;
}

sub dbID                   { return shift->_generic_get_or_set('dbID',                   @_) }
sub db                     { return shift->_generic_get_or_set('db',                     @_) }

=head2 name

  Example    : my $name = $regulatory_build->name;
  Description: Getter of name attribute for Regulatory Build objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name                   { return shift->_generic_get_or_set('name',                   @_) }

=head2 version

  Example    : my $name = $regulatory_build->version;
  Description: Getter of version attribute for Regulatory Build objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub version                { return shift->_generic_get_or_set('version',                @_) }

=head2 initial_release_date

  Example    : my $name = $regulatory_build->initial_release_date;
  Description: Getter of initial release date attribute for Regulatory Build objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub initial_release_date   { return shift->_generic_get_or_set('initial_release_date',   @_) }

=head2 last_annotation_update

  Example    : my $name = $regulatory_build->last_annotation_update;
  Description: Getter of last update attribute for Regulatory Build objects
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub last_annotation_update { return shift->_generic_get_or_set('last_annotation_update', @_) }
sub feature_type_id        { return shift->_generic_get_or_set('feature_type_id',        @_) }
sub analysis_id            { return shift->_generic_get_or_set('analysis_id',            @_) }

sub set_Analysis {
  my $self     = shift;
  my $analysis = shift;
  
  $self->{analysis} = $analysis;

  if (! defined $analysis) {
    throw('Analysis was not defined!');
  }
  
  $self->analysis_id($analysis->dbID);
}

sub fetch_Analysis {
  my $self = shift;
  
  if ($self->{analysis}) {
    return $self->{analysis};
  }
  my $analysis_id = $self->analysis_id;
  if ($analysis_id) {
    my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysis_id);
    $self->{analysis} = $analysis;
    return $analysis;
  }
  die;
}

sub set_FeatureType {
  my $self         = shift;
  my $feature_type = shift;
  
  if (! defined $feature_type) {
    throw('Feature type was not defined!');
  }
  $self->{feature_type} = $feature_type;
  $self->feature_type_id($feature_type->dbID);
}

sub fetch_FeatureType {
  my $self = shift;
  
  if ($self->{feature_type}) {
    return $self->{feature_type};
  }
  my $feature_type_id = $self->feature_type_id;
  if ($feature_type_id) {
    my $feature_type_adaptor = $self->db->get_FeatureTypeAdaptor;
    my $feature_type = $feature_type_adaptor->fetch_by_dbID($feature_type_id);
    $self->{feature_type} = $feature_type;
    return $feature_type;
  }
  die;
}

=head2 is_current

  Example    : if ($regulatory_build->is_current) { 
                 print "The current regulatory build was released on: " . $regulatory_build->initial_release_date . "\n";
               };
  Description: Indicates whether this regulatory build is the current one.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub is_current             { return shift->_generic_get_or_set('is_current',             @_) }

=head2 get_all_Epigenomes

  Example    : my $epigenome_in_regulatory_build = $regulatory_build->get_all_Epigenomes;
  Description: Gets all epigenomes used in this regulatory build.
  Returntype : ArrayRef[Bio::EnsEMBL::Funcgen::Epigenome]
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_all_Epigenomes {
  my $self = shift;
  
  my $epigenome_ids = $self->_get_all_epigenome_ids;
  my $epigenome_adaptor = $self->db->get_EpigenomeAdaptor;
  my $epigenomes = $epigenome_adaptor->fetch_all_by_dbID_list($epigenome_ids);
  return $epigenomes
}

sub _get_all_epigenome_ids {
  my $self = shift;
  
  my $db = $self->db;
  
  my $regulatory_build_id = $self->dbID;
  
  my $sql = qq(select epigenome_id from regulatory_build_epigenome where regulatory_build_id = ) . $regulatory_build_id;
  
  my $sth = $db->dbc->prepare( $sql );
  my $rv  = $sth->execute();
  my $res = $sth->fetchall_arrayref;
  
  my @epigenome_id = map { $_->[0] } @$res;

  return \@epigenome_id
}

=head2 fetch_sample_regulatory_feature

  Example    : print "Stable id of sample regulatory feature: " . $regulatory_build->fetch_sample_RegulatoryFeature->stable_id . "\n";
  Description: Gets all epigenomes used in this regulatory build.
  Returntype : ArrayRef[Bio::EnsEMBL::Funcgen::Epigenome]
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub fetch_sample_RegulatoryFeature {
  my $self = shift;
  
  my $sample_regulatory_feature_id = $self->sample_regulatory_feature_id;
  return if (! defined $sample_regulatory_feature_id);
  
  my $regulatory_feature_adaptor = $self->db->get_RegulatoryFeatureAdaptor;
  my $sample_regulatory_feature = $regulatory_feature_adaptor->fetch_by_dbID($sample_regulatory_feature_id);
  return $sample_regulatory_feature
}

sub sample_regulatory_feature_id { return shift->_generic_get_or_set('sample_regulatory_feature_id', @_) }

1;


