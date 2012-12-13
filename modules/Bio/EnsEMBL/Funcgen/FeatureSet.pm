#
# Ensembl module for Bio::EnsEMBL::Funcgen::FeatureSet
#

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::Funcgen::FeatureSet - A container for a discrete set of experimental features.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::FeatureSet;

my $feature_set = Bio::EnsEMBL::Funcgen::FeatureSet->new
  (
   -dbid          => $dbid,
   -analysis      => $analysis,
   -feature_type  => $ftype,
   -cell_type     => $ctype,
   -name          => 'HeLa-S3_H3K4me3_ExperimentalGroup_Analysis',
   -feature_class => 'annotated',
   -description   => 'HeLa-S3 H3K4me3 SWEmbl peaks replicate 1',
   -display_label => 'HeLa-S3 H3K4me3',
   -input_set     => $iset,
  ); 


my @features = @{$feature_set->get_Features_by_Slice($slice)};

print $feature_set->display_label.' has '.
  scalar(@features).' on slice '.$slice->name."\n";


$feature_set        = $feature_set_adaptor->fetch_by_name('RegulatoryFeatures:HepG2');
my $feature_type    = $feature_type_adaptor->fetch_by_name('Promoter Associated');

#This may take a while and a lot of memory.
my @prom_like_regfs = @{$feature_set->get_all_by_FeatureType($feature_type)};

#This may take longer and more memory.
my @all_regfs       = @{$feature_set->get_all_Features()};

=head1 DESCRIPTION

A FeatureSet object defines a discrete set of features from a single analysis.
It is a container object providing access to meta data and convenience wrapper 
methods to fetch the related features.

A FeatureSet is normally associated to a DataSet either as a product FeatureSet,
with the excpetion of FeatureSets that define external data. As these generally 
have no supporting data, there is no need for a DataSet. FeatureSets can also be 
associated with DataSets as supporting sets i.e. the input to a further analysis
such as the Ensembl Regulatory Build:
  http://www.ensembl.org/info/docs/funcgen/regulatory_build.html

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::Set;
Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor;
Bio::EnsEMBL::Funcgen::DataSet;

=cut

package Bio::EnsEMBL::Funcgen::FeatureSet;
use base qw(Bio::EnsEMBL::Funcgen::Set);
#use base rather than @ISA as we don't need to export from Set

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate);
use Bio::EnsEMBL::Funcgen::Set;

use strict;
use warnings;

my %valid_classes = (
                     annotated    => undef,
                     regulatory   => undef,
                     external     => undef,
                     segmentation => undef,
                    );

=head2 new

  Args [1] : Hash of key value parameter as follows:

             MANDATORY:
               -NAME          => String - name for this Set.
               -FEATURE_TYPE  => Bio::EnsEMBL::Funcgen::FeatureType
               -FEATURE_CLASS => String - Class of feature e.g. result, annotated,
                                 regulatory, segmentation, external or dna_methylation.
             OPTIONAL:
               -INPUT_SET     => Bio::EnsEMBL::Funcgen::InputSet
               -DESCRIPTION   => String
               -DISPLAY_NAME  => String
               -CELL_TYPE     => Bio::EnsEMBL::Funcgen::CellType
               -ANALYSIS      => Bio::EnsEMBL::Analysis
               -DBID          => Int
               -ADAPTOR       => Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor

  Example    : my $feature = Bio::EnsEMBL::Funcgen::FeatureSet->new
                 (
                  -dbid          => $dbid,
                  -analysis      => $analysis,
                  -feature_type  => $ftype,
                  -cell_type     => $ctype,
                  -name          => 'HeLa-S3_H3K4me3_ExperimentalGroup_Analysis',
                  -feature_class => 'annotated',
                  -description   => 'HeLa-S3 H3K4me3 SWEmbl peaks replicate 1',
                  -display_label => 'HeLa-S3 H3K4me3',
                  -input_set     => $iset,
			           );
  Description: Constructor for FeatureSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : Throws if: FeatureType is not defined, Feature class is not valid
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = $class->SUPER::new(@_);

  my ( $desc, $dlabel, $iset_id, $iset, $exp_id, $exp ) =
    rearrange( [
                'DESCRIPTION',   'DISPLAY_LABEL',
                'INPUT_SET_ID', 'INPUT_SET', 'EXPERIMENT_ID', 'EXPERIMENT'
               ], @_ );


  if ($exp_id || $exp) {
    throw('Passing an Experiment or an experiment_id is now deprecated,'.
          ' please use -input_set or -input_set_id instead');
  }

  #Allow exp or exp_id to be passed to support storing and lazy loading

  #Mandatory params checks here (setting done in Set.pm)
  if ( ! defined $self->feature_type ) {
    throw('Must provide a FeatureType');
  }

  #explicit type check here to avoid invalid types being imported as NULL
  #subsequently throwing errors on retrieval
  my $type = $self->feature_class;

  if ( ! ( $type && exists $valid_classes{$type} ) ) {
    throw( 'You must define a valid FeatureSet type e.g. ' .
           join( ', ', keys %valid_classes ) );
  }

  #Direct assignment to prevent need for set arg test in method

  $self->{description}   = $desc    if defined $desc;
  $self->{display_label} = $dlabel  if defined $dlabel;
  $self->{input_set_id}  = $iset_id if defined $iset_id;

  if ( defined $iset ) {
    #Exp obj is only passed during object storing
    #so let the adaptor do is_stored_and_valid
    $self->{input_set} = $iset;
  }

  return $self;
}                               ## end sub new



=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub new_fast { return bless ($_[1], $_[0]); }


=head2 description

  Example    : print "Feature set description is:\t".$fset->description."\n";
  Description: Getter for the description of this FeatureSet. e.g. Release 3.1
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub description { return $_[0]->{description}; }



=head2 display_label

  Example    : print $rset->display_label;
  Description: Getter for the display_label attribute for this FeatureSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
  my ($self) = @_;

  if ( ! defined $self->{display_label} ) {

    if ( $self->feature_class eq 'regulatory' ) {
      $self->{display_label} = $self->name;
    } else {
      #This still fails here if we don't have a class or a cell_type set

      $self->{display_label} =
        $self->feature_type->name.' - '. $self->cell_type->name.' enriched sites';
    }
  }

  return $self->{display_label};
}




=head2 get_FeatureAdaptor

  Example    :
  Description: Retrieves and caches FeatureAdaptor of feature_set type
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_FeatureAdaptor{
  my $self = shift;

  if (! exists $self->{'adaptor_refs'}) {

    foreach my $valid_class (keys %valid_classes) {
      my $method = 'get_'.$self->adaptor->build_feature_class_name($valid_class).'Adaptor';
      $self->{'adaptor_refs'}{$valid_class} = $self->adaptor->db->$method;
    }
  }

  return $self->{'adaptor_refs'}->{$self->feature_class()};
}



=head2 get_Features_by_Slice

  Example    : my @features = @{$FeatureSet->get_Features_by_Slice($slice)};
  Description: Retrieves all Features for this FeatureSet for a given Slice
  Returntype : ARRAYREF containing Features of the feature_set type i.e. Annotated, Regulatory or Supporting;
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Features_by_Slice{
  my ($self, $slice) = @_;

  return $self->get_FeatureAdaptor->fetch_all_by_Slice_FeatureSets($slice, [$self]);
}


=head2 get_Features_by_FeatureType

  Arg[0]     : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : my @features = @{$FeatureSet->get_Features_by_FeatureType($ftype)};
  Description: Retrieves all Features for this FeatureSet for a given FeatureType
               or associated FeatureType. This is mainly used by external FeatureSets
               which can sometimes have more than one associated FeatureType.
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_Features_by_FeatureType{
  my ($self, $type) = @_;

  return $self->get_FeatureAdaptor->fetch_all_by_FeatureType_FeatureSets($type, [$self]);
}


=head2 get_all_Features

  Example    : my @features = @{$FeatureSet->get_all_Features};
  Description: Retrieves all Features for this FeatureSet
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_all_Features{
  my $self = shift;

  return $self->get_FeatureAdaptor->fetch_all_by_FeatureSets([$self]);
}


#focus is old terminology
#tend to use 'core' now
#althoughh here we use fouc to determine a set which is actually used in the build
#where as core is at the feature type level
#i.e. you can have a feature set with a core feature type
#which is not in the build and is therefore not a focus set.

=head2 is_focus_set

  Args       : None
  Example    : if($fset->is_focus_set){ ... }
  Description: Returns true if FeatureSet is a focus set used in the RegulatoryBuild
  Returntype : Boolean
  Exceptions : Throws if meta entry not present
  Caller     : General
  Status     : At Risk

=cut

sub is_focus_set{
  my $self = shift;

  if (! defined $self->{focus_set}) {

    if (! defined $self->cell_type) {
      warn "FeatureSet without an associated CellType cannot be a focus set:\t".$self->name;
      $self->{focus_set} = 0;
    } else {
      $self->{focus_set} = $self->adaptor->fetch_focus_set_config_by_FeatureSet($self);
    }
  }

  return $self->{focus_set};
}


=head2 is_attribute_set

  Args       : None
  Example    : if($fset->is_attribute_set){ ... }
  Description: Returns true if FeatureSet is a supporting/attribute(focus or not) set used in the RegulatoryBuild
  Returntype : Boolean
  Exceptions : Throws if meta entry not present
  Caller     : General
  Status     : At Risk

=cut

sub is_attribute_set{
  my $self = shift;

  if (! defined $self->{attribute_set}) {

    if (! defined $self->cell_type) {
      warn "FeatureSet without an associated CellType cannot be a attribute set:\t".$self->name;
      $self->{attribute_set} = 0;
    } else {
      $self->{attribute_set} = $self->adaptor->fetch_attribute_set_config_by_FeatureSet($self);
    }
  }

  return $self->{attribute_set};
}


=head2 get_InputSet

  Example    : my $input_set = $FeatureSet->get_InputSet;
  Description: Retrieves the InputSet for this FeatureSet
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_InputSet{
  my $self = shift;

  if ( (! defined $self->{input_set}) &&
       (defined $self->{input_set_id}) ) {
    $self->{input_set} = $self->adaptor->db->get_InputSetAdaptor->fetch_by_dbID($self->{input_set_id});
  }

  return $self->{input_set};
}



=head2 source_label

  Example    : my $source_label = $fset->source_label;
  Description: Retrieves the source label this FeatureSet, used in zmenus
  Returntype : Arrayref of Strings
  Exceptions : None
  Caller     : Webcode zmenus
  Status     : At Risk - remove, to be done by webcode? or move to InputSet and wrap from here

=cut

#These are used to link through to the experiment view based on feature_set_id
#select input_set_id, count(distinct archive_id) as cnt , group_concat(archive_id) from input_subset where archive_id is not NULL group by input_set_id having cnt >1;

sub source_label{
  my $self = shift;

  if (! defined $self->{source_label}) {
    my $input_set = $self->get_InputSet;
    my @source_labels;

    if ($input_set) {

      foreach my $isset (@{$input_set->get_InputSubsets}) {

        if (defined $isset->archive_id) {
          push @source_labels, $isset->archive_id;
        }
        #Archive IDs e.g. SRX identifiers or undef.
      }

      #Append project name
      my $exp_group = $input_set->get_Experiment->experimental_group;

      if ($exp_group &&
          $exp_group->is_project) {
        push @source_labels, $exp_group->name;
      }
    }

    $self->{source_label} = join(q{ }, # Single space
                                 @source_labels);
  }

  return $self->{source_label};
}




### DEPRECATED ###

=head2 get_Experiment
#
#  Example    : my $exp = $FeatureSet->get_Experiment;
#  Description: Retrieves the Experiment for this FeatureSet
#  Returntype : Bio::EnsEMBL::Funcgen::Experiment
#  Exceptions : None
#  Caller     : General
  Status     : DEPRECATED

=cut

sub get_Experiment{
  throw('FeatureSet::get_Experiment is not longer supported, please use FeatureSet::get_InputSet');
}


1;
