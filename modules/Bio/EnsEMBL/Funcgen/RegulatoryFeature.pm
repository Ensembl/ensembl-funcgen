=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

  Bio::EnsEMBL::Funcgen::RegulatoryFeature

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
  );

  my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryFeature');
  my $regulatory_feature = $regulatory_feature_adaptor->fetch_by_stable_id('ENSR00000000011');

  print 'Stable id:        ' . $regulatory_feature->stable_id                              . "\n";
  print 'Analysis:         ' . $regulatory_feature->analysis->logic_name                   . "\n";
  print 'Feature type:     ' . $regulatory_feature->feature_type->name                     . "\n";
  print 'Epigenome count:  ' . $regulatory_feature->epigenome_count                        . "\n";
  print 'Slice name:       ' . $regulatory_feature->slice->name                            . "\n";
  print 'Coordinates:      ' . $regulatory_feature->start .' - '. $regulatory_feature->end . "\n";
  print 'Regulatory build: ' . $regulatory_feature->get_regulatory_build->name             . "\n";

=head1 DESCRIPTION

A RegulatoryFeature object represents the output of the Ensembl RegulatoryBuild:
    http://www.ensembl.org/info/docs/funcgen/regulatory_build.html

It may comprise many histone modification, transcription factor, polymerase and open
chromatin features, which have been combined to provide a summary view and
classification of the regulatory status at a given loci.


=head1 SEE ALSO

Bio::EnsEMBL:Funcgen::DBSQL::RegulatoryFeatureAdaptor

=cut


package Bio::EnsEMBL::Funcgen::RegulatoryFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use base qw( Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  Arg [-SLICE]             : Bio::EnsEMBL::Slice - The slice on which this feature is located.
  Arg [-START]             : int - The start coordinate of this feature relative to the start of the slice
                             it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]               : int -The end coordinate of this feature relative to the start of the slice
                    	     it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-FEATURE_SET]       : Bio::EnsEMBL::Funcgen::FeatureSet - Regulatory Feature set
  Arg [-FEATURE_TYPE]      : Bio::EnsEMBL::Funcgen::FeatureType - Regulatory Feature sub type
  Arg [-STABLE_ID]         : (optional) string - Stable ID for this RegulatoryFeature e.g. ENSR00000000001
  Arg [-DISPLAY_LABEL]     : (optional) string - Display label for this feature
  Arg [-PROJECTED]         : (optional) boolean - Flag to specify whether this feature has been projected or not
  Arg [-dbID]              : (optional) int - Internal database ID.
  Arg [-ADAPTOR]           : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.

  Example    : my $feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new(
                -SLICE         => $chr_1_slice,
                -START         => 1000000,
                -END           => 1000024,
                -DISPLAY_LABEL => $text,
                -FEATURE_SET   => $fset,
                -FEATURE_TYPE  => $reg_ftype,
              );

  Description: Constructor for RegulatoryFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($stable_id, $attr_cache, $projected, $activity, $epigenome_count, $analysis)
    = rearrange(['STABLE_ID', 'PROJECTED', 'ACTIVITY', 'EPIGENOME_COUNT', 'ANALYSIS'], @_);

  #None of these are mandatory at creation
  #under different use cases
  $self->{stable_id}        = $stable_id        if defined $stable_id;
  $self->{projected}        = $projected        if defined $projected;
  $self->{activity}         = $activity         if defined $activity;
  $self->{epigenome_count}  = $epigenome_count  if defined $epigenome_count;
  $self->{analysis}         = $analysis         if defined $analysis;
  
  $self->{_regulatory_activity} = [];

  return $self;
}

sub analysis {
  return shift->get_Analysis;
}

sub _analysis_id {
  return shift->{_analysis_id};
}

=head2 get_Analysis

  Arg : none
  Example    : $analysis = $regulatory_feature->get_Analysis
  Description: Fetches the analysis used to generate this regulatory feature.
               This is the analysis of the regulatory build.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Status     : At Risk

=cut

sub get_Analysis {
  my $self = shift;

  if(! defined $self->{'_analysis'}) {
    $self->{'_analysis'} = $self
      ->adaptor
      ->db
      ->get_AnalysisAdaptor()
      ->fetch_by_dbID(
	$self->_analysis_id
      );
  }
  return $self->{'_analysis'};
}

sub regulatory_build_id {
  return shift->{regulatory_build_id};
}

=head2 display_label

  Example    : my $label = $feature->display_label;
  Description: Getter for the display label of this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_label {
  my $self = shift;

  if(! defined $self->{display_label}) {
    $self->{'display_label'}  = $self->feature_type->name.' Regulatory Feature';
  }

  return $self->{display_label};
}

sub feature_type {
  my $self = shift;
  my $value = shift;

  if(defined $value) {
    $self->{'feature_type'}  = $value;
  }

  return $self->{feature_type};
}

=head2 display_id

  Example    : print $feature->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. In this case it is the stable id of 
               the regulatory feature.
  Returntype : String
  Exceptions : none
  Caller     : web drawing code, Region Report tool
  Status     : Stable

=cut

sub display_id {  return shift->{stable_id}; }

=head2 stable_id

  Example    : my $stable_id = $feature->stable_id();
  Description: Getter for the stable_id attribute for this feature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At risk - setter functionality to be removed

=cut

sub stable_id { return shift->{stable_id}; }

=head2 regulatory_evidence

  Arg [1]    : String (optional) - Class of feature e.g. 'annotated' 
               or 'motif'
  Arg [2]    : Bio::EnsEMBL::Funcgen::Epigenome - The epigenome for which 
               the evidence for it regulatory activity is requested.
  Example    : 
  Description: Deprecated: Use get_RegulatoryEvidence instead.
  Returntype : ARRAYREF
  Exceptions : Throws if feature class not valid
  Caller     : General
  Status     : Deprecated

=cut

sub regulatory_evidence {
  deprecate(
    "Bio::EnsEMBL::Funcgen::RegulatoryFeature::regulatory_evidence() has been deprecated and will be removed in Ensembl release 91."
        . " Please use Bio::EnsEMBL::Funcgen::RegulatoryFeature::get_RegulatoryEvidence() instead."
);

  my $self = shift;
  my $feature_class = shift;
  my $epigenome   = shift;
  
  return $self->get_RegulatoryEvidence($feature_class, $epigenome);
}

=head2 get_RegulatoryEvidence

  Arg [1]    : String - Class of feature e.g. 'annotated'
               or 'motif'
  Arg [2]    : Bio::EnsEMBL::Funcgen::Epigenome - The epigenome for which 
               the evidence for its regulatory activity is requested.
               
               If no epigenome is provided, it will return a summary that
               can be used for the regulatory build track on the Ensembl
               website.
  Example    : 
  Description: Getter for the regulatory_evidence for this feature.
  Returntype : ARRAYREF
  Exceptions : Throws if feature class not valid
  Caller     : General
  Status     : At Risk

=cut

sub get_RegulatoryEvidence {
  my $self = shift;
  my $feature_class = shift;
  my $epigenome   = shift;

  if(! defined $feature_class){
    throw("Feature class string ('annotated' or 'motif') not defined!");
  }
  
  if (! defined $epigenome) {
	# Hack, so we get a summary for the regulatory build track until the
	# current data issues have been fixed.
	$epigenome = $self->adaptor->db->get_EpigenomeAdaptor->fetch_by_dbID(1);
  }
  $self->_assert_epigenome_ok($epigenome);
  my $regulatory_activity = $self->regulatory_activity_for_epigenome($epigenome);
  
  # See https://github.com/Ensembl/ensembl-funcgen/pull/6
  return [] unless $regulatory_activity;
  
  my $regulatory_evidence = $regulatory_activity->get_RegulatoryEvidence;

  return $regulatory_evidence;
}

sub _assert_epigenome_ok {
  my $self = shift;
  my $epigenome = shift;
  if (! defined $epigenome) {
    throw("Epigenome parameter was undefined!");
  }
  if (ref $epigenome ne 'Bio::EnsEMBL::Funcgen::Epigenome') {
    throw("epigenome parameter must have type Bio::EnsEMBL::Funcgen::Epigenome!");
  }
}

=head2 regulatory_activity_for_epigenome

  Arg [1]    : Bio::EnsEMBL::Funcgen::Epigenome - The epigenome for which 
               the evidence for it regulatory activity is requested.
  Example    : 
  Description: Getter for the regulatory_activity_for_epigenome for this 
               feature. Returns undef, if no activity was predicted for the
               epigenome. 
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryActivity
  Exceptions : none
  Caller     : General
  Status     : At Risk

=cut

sub regulatory_activity_for_epigenome {
  my $self = shift;
  my $epigenome = shift;

  if (! defined $epigenome) {
    throw("Epigenome parameter was undefined!");
  }
  if (ref $epigenome ne 'Bio::EnsEMBL::Funcgen::Epigenome') {
    throw("Wrong parameter, expected an epigenome, but got a " . ref $epigenome);
  }
  
  my $epigenome_id = $epigenome->dbID;
  my @regulatory_activity = grep { 
    !$_->_is_multicell 
    && $_->get_Epigenome->dbID == $epigenome_id 
  } @{$self->regulatory_activity};
  
  if (! @regulatory_activity) {
    return;
  }
  if (@regulatory_activity>1) {
    throw();
  }
  return $regulatory_activity[0];
}

sub _get_underlying_structure_motifs_by_epigenome {
  my $self = shift;
  my $epigenome = shift;
  
  my $regulatory_activity = $self->regulatory_activity_for_epigenome($epigenome);
  
  if (! defined $regulatory_activity) {
    die("No regulatory activity found for epigenome " . $epigenome->display_label . ".");
  }
  my $motif_evidence = $regulatory_activity->get_RegulatoryEvidence_by_type('motif');
  
  return $motif_evidence;
}

sub _get_underlying_structure_motifs_by_epigenome_list {
  my $self = shift;
  my $epigenome_list = shift;
  
  my @all_motifs;
  
  foreach my $current_epigenome (@$epigenome_list) {
    my $motif_evidence = $self->_get_underlying_structure_motifs_by_epigenome($current_epigenome);
    push @all_motifs, @$motif_evidence;
  }
  
  my @unique_motifs;
  my %seen_motif_coordinates;

  MOTIF: foreach my $current_motif (@all_motifs) {
  
    my $unique_key = $current_motif->start . '_' . $current_motif->end;
    if (exists $seen_motif_coordinates{$unique_key}) {
      next MOTIF;
    }
    push @unique_motifs, $current_motif;
    $seen_motif_coordinates{$unique_key} = 1;
  }
  my @sorted_motifs = sort { $a->start <=> $b->start } @unique_motifs;
  return \@sorted_motifs;
}

=head2 get_underlying_structure

  Example    : my @web_image_structure = @{$regulatory_feature->get_underlying_structure};
  
  Description: Returns the underlying structure of the regulatory feature in a
               given epigenome. 

               This structure is the bound start and end and the boundaries 
               of all motifs linked to this regulatory feature for this 
               epigenome. 

               If no epigenome is provided, it will return a summary, which 
               is the union of all underlying structures of all epigenomes in 
               the regulatory build.

  Returntype : Arrayref of slice coordinates relative to the slice of the 
               regulatory feature.
  Exceptions : None
  Caller     : Webcode
  Status     : At Risk

=cut

sub get_underlying_structure {
  my $self = shift;
  my $epigenome = shift;
  
  my $epigenome_specific_underlying_structure = [];
  my $motif_evidence;
  
  if ($epigenome) {
    $self->_assert_epigenome_ok($epigenome);
    $motif_evidence = $self->_get_underlying_structure_motifs_by_epigenome($epigenome);
  } else {
    my $regulatory_build = $self->get_regulatory_build;
    my $epigenome_list = $regulatory_build->get_all_Epigenomes;
    $motif_evidence = $self->_get_underlying_structure_motifs_by_epigenome_list($epigenome_list);
  }
  
  # This is used to make the coordinates relative to the start of the current
  # slice. The webcode expects it that way.
  #
  my $slice_start = $self->slice->start;
  
  foreach my $current_motif_evidence (@$motif_evidence) {

    push @$epigenome_specific_underlying_structure, (
      0 + $current_motif_evidence->start - $slice_start +1, 
      0 + $current_motif_evidence->end   - $slice_start +1, 
    );
  }
  my $underlying_structure = [
    0 + $self->bound_start, 
    0 + $self->start,
    @$epigenome_specific_underlying_structure,
    0 + $self->end, 
    0 + $self->bound_end
  ];
  return $underlying_structure;
}

=head2 regulatory_activity

  Example    : my $regulatory_activity = $regulatory_feature->regulatory_activity;
  
  Description: Fetches a list of regulatory activities associated to this 
               regulatory feature.

  Returntype : ArrayRef[Bio::EnsEMBL::Funcgen::RegulatoryActivity]
  Exceptions : None
  Status     : At Risk

=cut

sub regulatory_activity {

  my $self = shift;
  if(! defined $self->{'_regulatory_activity'}) {
    my $raa = $self->adaptor->db->{'_regulatory_activity_adaptor'} ||= $self
      ->adaptor
      ->db
      ->get_RegulatoryActivityAdaptor();

    $self->{'_regulatory_activity'} = $raa->fetch_all_by_RegulatoryFeature($self);
  }
  return $self->{'_regulatory_activity'};
}

=head2 get_regulatory_build

  Example    : my $regulatory_build = $regulatory_feature->get_regulatory_build;
  
  Description: Fetches the regulatory build used to generate this regulatory 
               feature.

  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryBuild
  Exceptions : None
  Status     : At Risk

=cut

sub get_regulatory_build {
  my $self = shift;

  if(! defined $self->{'_regulatory_build'}) {
    $self->{'_regulatory_build'} = $self
      ->adaptor
      ->db
      ->get_RegulatoryBuildAdaptor()
      ->fetch_by_dbID(
	$self->regulatory_build_id
      );
  }
  return $self->{'_regulatory_build'};
}

sub add_regulatory_activity {

  my $self = shift;
  my $regulatory_activity = shift;
  
  push @{$self->{_regulatory_activity}}, $regulatory_activity;
}

sub has_activity_in {

  my $self = shift;
  my $epigenome = shift;
  
  foreach my $current_regulatory_activity (@{$self->regulatory_activity}) {
    if ($current_regulatory_activity->_epigenome_id() == $epigenome->dbID) {
      return 1;
    }
  }
  return;
}

sub has_epigenomes_with_activity {

  my $self = shift;
  my $activity = shift;
  
  foreach my $current_regulatory_activity (@{$self->regulatory_activity}) {
    if ($current_regulatory_activity->activity eq $activity) {
      return 1;
    }
  }
  return;
}

=head2 get_epigenomes_by_activity

  Example    : my $regulatory_build = $regulatory_feature->get_epigenomes_by_activity($regulatory_activity);
  
  Description: Returns an array of epigenomes in which this regulatory 
               feature has the regulatory activity given to this method as
               parameter.
               
               This can be used to get a list of epigenomes in which a 
               regulatory feature is active.

  Returntype : ArrayRef[Bio::EnsEMBL::Funcgen::Epigenome]
  Exceptions : None
  Status     : At Risk

=cut

sub get_epigenomes_by_activity {

  my $self     = shift;
  my $activity = shift;
  
  if (!$self->adaptor->is_valid_activity($activity)) {
    throw(
      'Please pass a valid activity to this method. Valid activities are: ' 
      . $self->adaptor->valid_activities_as_string
    );
  }
  
  my @epigenome_dbID_list = grep {
  
    # Multicell does not have an epigenome id. In the map statement below 
    #
    #   $_->epigenome_id 
    #
    # will come up with undef and the undef is removed here.
    
    $_ 
  } map { 
    $_->_epigenome_id 
  } grep { 
    $_->activity eq $activity 
  } @{$self->regulatory_activity};
  
  my $epigenome_adaptor = $self->adaptor->db->get_EpigenomeAdaptor;
  
  return $epigenome_adaptor->fetch_all_by_dbID_list(\@epigenome_dbID_list);
}

=head2 epigenome_count

  Arg [1]     : None
  Returntype  : SCALAR
  Exceptions  : None
  Description : Returns the amount of epigenomes in which this regulatory feature is active

=cut

sub epigenome_count { 
  my $self = shift;
  return $self->{epigenome_count};
}

=head2 bound_seq_region_start

  Example    : my $bound_sr_start = $feature->bound_seq_region_start;
  Description: Getter for the seq_region bound_start attribute for this feature.
               Gives the 5' most start value of the underlying attribute
               features.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_seq_region_start { return $_[0]->seq_region_start - $_[0]->_bound_lengths->[0]; }

=head2 bound_seq_region_end

  Example    : my $bound_sr_end = $feature->bound_seq_region_end;
  Description: Getter for the seq_region bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_seq_region_end { return $_[0]->seq_region_end + $_[0]->_bound_lengths->[1]; }

sub _bound_lengths {
  my $self = shift;

  if(! defined  $self->{_bound_lengths}){

    my @af_attrs = @{$self->regulatory_evidence('annotated')};

    if (! @af_attrs) {
      throw('Unable to set bound length, no AnnotatedFeature attributes available for RegulatoryFeature: '
            .$self->dbID);
    }

    #Adding self here accounts for core region i.e.
    #features extending beyond the core may be absent on this cell type.
    my @start_ends;

    foreach my $feat (@af_attrs, $self) {
      push @start_ends, ($feat->seq_region_start, $feat->seq_region_end);
    }

    @start_ends = sort { $a <=> $b } @start_ends;

    $self->{_bound_lengths} = [ ($self->seq_region_start - $start_ends[0]),
                                ($start_ends[$#start_ends] - $self->seq_region_end) ];
  }

  return $self->{_bound_lengths};
}

=head2 bound_start_length

  Example    : my $bound_start_length = $reg_feat->bound_start_length;
  Description: Getter for the bound_start_length attribute for this feature,
               with respect to the host slice strand
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_start_length {
  my $self = shift;
  return ($self->slice->strand == 1) ? $self->_bound_lengths->[0] : $self->_bound_lengths->[1];
}


=head2 bound_end_length

  Example    : my $bound_end_length = $reg_feat->bound_end_length;
  Description: Getter for the bound_end length attribute for this feature,
               with respect to the host slice strand.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub bound_end_length {
  my $self = shift;
  return ($self->slice->strand == 1) ? $self->_bound_lengths->[1] : $self->_bound_lengths->[0];
}

=head2 bound_start

  Example    : my $bound_start = $feature->bound_start;
  Description: Getter for the bound_start attribute for this feature.
               Gives the 5' most start value of the underlying attribute
               features in local coordinates.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub bound_start { return $_[0]->start - $_[0]->bound_start_length; }


=head2 bound_end

  Example    : my $bound_end = $feature->bound_start();
  Description: Getter for the bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features in local coordinates.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub bound_end { return $_[0]->end + $_[0]->bound_end_length; }

=head2 is_projected

  Arg [1]    : optional - boolean
  Example    : if($regf->is_projected){ #do something different here }
  Description: Getter/Setter for the projected attribute.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : At risk - remove setter functionality

=cut
sub is_projected {
  my $self = shift;
  deprecate("Bio::EnsEMBL::Funcgen::RegulatoryFeature::is_projected() has " .
  "been deprecated and will be removed in Ensembl release 93");

  if(@_){
	#added v67
    warn "RegulatoryFeature::is_projected setter functionality is being removed\n";
    $self->{'projected'} = shift;
  }

  return $self->{'projected'};
}

=head2 summary_as_hash

  Example       : $regulatory_feature_summary = $regulatory_feature->summary_as_hash;
  Description   : Retrieves a textual summary of this RegulatoryFeature.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self   = shift;
  
  my $feature_type = $self->feature_type;

  return {
    ID                => $self->stable_id,
    source            => $self->analysis->logic_name,,
    bound_start       => $self->bound_seq_region_start,
    bound_end         => $self->bound_seq_region_end,
    start             => $self->seq_region_start,
    end               => $self->seq_region_end,
    strand            => $self->strand,
    seq_region_name   => $self->seq_region_name,
    description       => $self->feature_type->description,
    feature_type      => $feature_type->name,
  };
}

1;


