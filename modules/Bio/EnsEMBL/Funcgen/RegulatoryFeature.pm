=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

  use v5.10;
  use Bio::EnsEMBL::Registry;
  Bio::EnsEMBL::Registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org:3306/84/');

  my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor('human', 'funcgen', 'RegulatoryFeature');

  my $regulatory_feat = $regulatory_feature_adaptor->fetch_by_stable_id('ENSR00000000021');

  say 'Stable id:       ' . $regulatory_feat->stable_id;
  say 'Analysis:        ' . $regulatory_feat->analysis->logic_name;
  say 'Feature type:    ' . $regulatory_feat->feature_type->name;
  say 'Feature set:     ' . $regulatory_feat->epigenome->name;
  say 'Activity:        ' . $regulatory_feat->activity;
  say 'Cell type count: ' . $regulatory_feat->cell_type_count;
  say 'Slice name:      ' . $regulatory_feat->slice->name;
  say 'Coordinates:     ' . $regulatory_feat->start .' - '. $regulatory_feat->end,;

=head1 DESCRIPTION

A RegulatoryFeature object represents the output of the Ensembl RegulatoryBuild:
    http://www.ensembl.org/info/docs/funcgen/regulatory_build.html

It may comprise many histone modification, transcription factor, polymerase and open
chromatin features, which have been combined to provide a summary view and
classification of the regulatory status at a given loci.


=head1 SEE ALSO

Bio::EnsEMBL:Funcgen::DBSQL::RegulatoryFeatureAdaptor
Bio::EnsEMBL::Funcgen::SetFeature

=cut


package Bio::EnsEMBL::Funcgen::RegulatoryFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use base qw( Bio::EnsEMBL::Funcgen::SetFeature );

=head2 new

  Arg [-SLICE]             : Bio::EnsEMBL::Slice - The slice on which this feature is located.
  Arg [-START]             : int - The start coordinate of this feature relative to the start of the slice
                             it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]               : int -The end coordinate of this feature relative to the start of the slice
                    	     it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-FEATURE_SET]       : Bio::EnsEMBL::Funcgen::FeatureSet - Regulatory Feature set
  Arg [-FEATURE_TYPE]      : Bio::EnsEMBL::Funcgen::FeatureType - Regulatory Feature sub type
  Arg [-BINARY_STRING]     : (optional) string - Regulatory Build binary string
  Arg [-STABLE_ID]         : (optional) string - Stable ID for this RegulatoryFeature e.g. ENSR00000000001
  Arg [-DISPLAY_LABEL]     : (optional) string - Display label for this feature
  Arg [-ATTRIBUTE_CACHE]   : (optional) HASHREF of feature class dbID|Object lists
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

  my ($stable_id, $attr_cache, $bin_string, $projected, $activity, $epigenome_count, $analysis)
    = rearrange(['STABLE_ID', 'ATTRIBUTE_CACHE', 'BINARY_STRING', 'PROJECTED', 'ACTIVITY', 'EPIGENOME_COUNT', 'ANALYSIS'], @_);

  #None of these are mandatory at creation
  #under different use cases
  $self->{binary_string}    = $bin_string       if defined $bin_string;
  $self->{stable_id}        = $stable_id        if defined $stable_id;
  $self->{projected}        = $projected        if defined $projected;
  $self->{activity}         = $activity         if defined $activity;
  $self->{epigenome_count}  = $epigenome_count  if defined $epigenome_count;
  $self->{epigenome_count}  = $epigenome_count  if defined $epigenome_count;
  $self->{analysis}         = $analysis         if defined $analysis;
  
  $self->{_regulatory_activity} = [];

  return $self;
}

# deprecated, use get_Analysis
sub analysis {
  return shift->get_Analysis;
}

sub _analysis_id {
  return shift->{_analysis_id};
}

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

    if( defined $self->epigenome ) {
      $self->{display_label} .= ' - '.$self->epigenome->name;
    }
  }

  return $self->{display_label};
}


=head2 display_id

  Example    : print $feature->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. In this case the stable Id is
               preferred
  Returntype : String
  Exceptions : none
  Caller     : web drawing code, Region Report tool
  Status     : Stable

=cut

sub display_id {  return shift->{stable_id}; }

=head2 stable_id

  Arg [1]    : (optional) string - stable_id e.g ENSR00000000001
  Example    : my $stable_id = $feature->stable_id();
  Description: Getter for the stable_id attribute for this feature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At risk - setter functionality to be removed

=cut

sub stable_id { return shift->{stable_id}; }


=head2 regulatory_evidence

  Arg [1]    : String (optional) - Class of feature e.g. annotated or motif
  Example    : print "Regulatory Attributes:\n\t".join("\n\t", (map $_->feature_type->name, @{$feature->regulatory_evidence()}))."\n";
  Description: Getter for the regulatory_evidence for this feature.
  Returntype : ARRAYREF
  Exceptions : Throws if feature class not valid
  Caller     : General
  Status     : At Risk

=cut
sub regulatory_evidence {
  my $self = shift;
  my $feature_class = shift;
  my $epigenome   = shift;
  
  $self->_assert_epigenome_ok($epigenome);
  my $regulatory_activity = $self->regulatory_activity_for_epigenome($epigenome);
  
  # See https://github.com/Ensembl/ensembl-funcgen/pull/6
  return [] unless $regulatory_activity;
  
  my $regulatory_evidence = $regulatory_activity->regulatory_evidence;
    
  if ($feature_class eq 'annotated') {
    return $regulatory_evidence->supporting_annotated_features;
  }
  if ($feature_class eq 'motif') {
    return $regulatory_evidence->supporting_motif_features;
  }
  throw("Invalid feature class $feature_class!");
}

sub _assert_epigenome_ok {
  my $self = shift;
  my $epigenome = shift;
  if (! defined $epigenome) {
    throw();
  }
  if (ref $epigenome ne 'Bio::EnsEMBL::Funcgen::Epigenome') {
    throw("epigenome parameter must have type Bio::EnsEMBL::Funcgen::Epigenome!");
  }
}

sub regulatory_activity_for_epigenome {
  my $self = shift;
  my $epigenome = shift;
  
  my $epigenome_id = $epigenome->dbID;
  my @regulatory_activity = grep { 
    !$_->_is_multicell 
    && $_->epigenome_id == $epigenome_id 
  } @{$self->regulatory_activity};
  
  if (! @regulatory_activity) {
#     throw('No regulatory activity for ' . $epigenome->display_label);
    return;
  }
  if (@regulatory_activity>1) {
    throw();
  }
  return $regulatory_activity[0];
}

=head2 get_underlying_structure

  Example    : my @web_image_structure = @{$regf->get_underlying_structure};
  Description: Getter for the bound_end attribute for this feature.
               Gives the 3' most end value of the underlying attribute
               features.
  Returntype : Arrayref
  Exceptions : None
  Caller     : Webcode
  Status     : At Risk

=cut

sub get_underlying_structure {
  my $self = shift;
  my $epigenome = shift;
#   $self->_assert_epigenome_ok($epigenome);
  

  # Stopgap fix to prevent error messages for the missing 
  # get_underlying_structure in regulatory_activity.
  #
  # Must have been lost in a git merge black hole.
  #
  # Hopefully resurrected soon and the real code below can be put in place 
  # again.
  #
  my $underlying_structure = [
    0 + $self->bound_start, 
    0 + $self->start,
    0 + $self->end, 
    0 + $self->bound_end
  ];

#   my $regulatory_activity = $self->regulatory_activity_for_epigenome($epigenome);
# 
#   my $epigenome_specific_underlying_structure = $regulatory_activity->get_underlying_structure();
# 
#   my $underlying_structure = [
#     0 + $self->bound_start, 
#     0 + $self->start,
#     @$epigenome_specific_underlying_structure,
#     0 + $self->end, 
#     0 + $self->bound_end
#   ];
# 
  return $underlying_structure;
}

=head2 regulatory_activity

  Arg [1]     : 
  Returntype  : 
  Exceptions  : 
  Description : Guaranteed to return an arrayref. If there are no linked feature sets, returns [].

=cut
sub regulatory_activity {

  my $self = shift;
  if(! defined $self->{'_regulatory_activity'}) {
    $self->{'_regulatory_activity'} = $self
      ->adaptor
      ->db
      ->get_RegulatoryActivityAdaptor()
      ->fetch_all_by_RegulatoryFeature($self);
  }
  return $self->{'_regulatory_activity'};
}

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
    if ($current_regulatory_activity->epigenome_id == $epigenome->dbID) {
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

  Arg [1]     : Activity
  Returntype  : 
  Exceptions  : 
  Description : 

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
    $_->epigenome_id 
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

  if(@_){
	#added v67
    warn "RegulatoryFeature::is_projected setter functionality is being removed\n";
    $self->{'projected'} = shift;
  }

  return $self->{'projected'};
}

=head2 summary_as_hash

  Example       : $regf_summary = $regf->summary_as_hash;
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

# Deprecated methods

sub has_evidence {
    deprecate('"has_evidence" is now deprecated. Please use "activity"
        which reports the state of the Regulatory Feature');
  return shift->activity;
}

sub cell_type_count { 
  my $self = shift;
  deprecate(
        "Bio::EnsEMBL::Funcgen::RegulatoryFeature::cell_type_count has been deprecated and will be removed in Ensembl release 89."
            . " Please use Bio::EnsEMBL::Funcgen::RegulatoryFeature::epigenome_count instead"
  );
  return $self->epigenome_count;
}

sub is_unique_to_FeatureSets { deprecate('"is_unique_to_FeatureSets" is deprecated. '); return; }
sub get_other_RegulatoryFeatures { deprecate('"get_other_RegulatoryFeatures" is deprecated. '); return; }
sub get_focus_attributes    { deprecate('"get_focus_attributes" is deprecated.');  return; }
sub get_nonfocus_attributes { deprecate('"get_nonfocus_attributes" is deprecated.');  return; }

sub activity {
  throw(
    "activity is no longer supported for regulatory features. You can use "
    . "get_epigenomes_by_activity('ACTIVE') to find feature sets in which "
    . "this regulatory feature are active."
  );
}

sub epigenome {
  throw(
    "epigenome is no longer supported for regulatory features. You can use "
    . "get_epigenomes_by_activity('ACTIVE') to find feature sets in which "
    . "this regulatory feature are active."
  );
}

1;


