#
# Ensembl module for Bio::EnsEMBL::Funcgen::SetFeature
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

Bio::EnsEMBL::Funcgen::SetFeature - Ensembl specific set feature.

=head1 SYNOPSIS

    my $feat = new Bio::EnsEMBL::Feature(-start         => 100,
                                         -end           => 220,
                                         -strand        => -1,
                                         -slice         => $slice,
                                         -feature_set   => $fset,
                                         -display_label => $label,
                                      );

    my $start  = $feat->start;
    my $end    = $feat->end;
    my $strand = $feat->strand;

    #move the feature to the chromosomal coordinate system
    $feature = $feature->transform('chromosome');

    #move the feature to a different slice (possibly on another coord system)
    $feature = $feature->transfer($new_slice);

    #project the feature onto another coordinate system possibly across
    #boundaries:
    @projection = @{$feature->project('contig')};

    #change the start, end, and strand of the feature in place
    $feature->move($new_start, $new_end, $new_strand);

=head1 DESCRIPTION

This is a simple wrapper method for the core Feature class to contain common generic
Funcgen SetFeature methods.

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::SetFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Funcgen::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable);


#Bio::EnsEMBL::Funcgen::Storable
#gives access to DBEntry type methods
#This only works as the DBEntry methods in the core API are at the Gene/Transcript/Translation level
#This would not be a problem if this was integrated back into the core code

=head2 new

  Arg [-SET]          : Bio::EnsEMBL::Funcgen::ResultSet or FeatureSet.
  Arg [-FEATURE_TYPE] : optional Bio::EnsEMBL::Funcgen::FeatureType. Defaults to Set FeatureType.
  Arg [-ANALYSIS]     : Bio::EnsEMBL::Analysis. Defaults to Set Analysis.
  Arg [-SLICE]        : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-START]        : int - The start coordinate of this feature relative to the start of the slice
		                it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int -The end coordinate of this feature relative to the start of the slice
	                    it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-DISPLAY_LABEL]: string - Display label for this feature
  Arg [-STRAND]       : int - The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]         : (optional) int - Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.

  Arg [-FEATURE_SET]  : Bio::EnsEMBL::Funcgen::FeatureSet (Obsolete)

  Example             : my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
                              (
                               -SLICE         => $chr_1_slice,
                               -START         => 1000000,
                               -END           => 1000024,
                               -STRAND        => -1,
							   -DISPLAY_LABEL => $text,
							   -SET           => $fset,
                              );

  Description: Constructor for SetFeature objects. Should never be called directly, only by children.
  Returntype : Bio::EnsEMBL::Funcgen::SetFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk - FEATURE_SET arg replaced by SET in v67

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;

  my ($display_label, $fset, $ftype, $set)
    = rearrange(['DISPLAY_LABEL', 'FEATURE_SET', 'FEATURE_TYPE', 'SET'], @_);
  
  #Don't test for set vs fset here as this will slow things down.
  #We use new_fast from the adaptors to avoid this
  $set ||= $fset;

  #Grab FeatureSet first so we can pass analysis to base Feature class
  if( ( ref($set) ne 'Bio::EnsEMBL::Funcgen::FeatureSet') &&
	  ( ref($set) ne 'Bio::EnsEMBL::Funcgen::ResultSet') ){
	throw("Must pass valid Bio::EnsEMBL::Funcgen::FeatureSet or ResultSet object");
  }

  my $self = $class->SUPER::new(@_, -analysis => $fset->analysis);
  
   

  if($ftype){
	
	if(ref($ftype) ne 'Bio::EnsEMBL::Funcgen::FeatureType'){
	  throw('feature_type param must be a valid Bio::EnsEMBL::Funcgen::FeatureType');
	}
  
	$self->{'feature_type'} = $ftype;
  }
 
  #Setting attrs directly removes the need for setter code in methods
  $self->{set}           = $set;
  $self->{display_label} = $display_label if defined $display_label;
 	
  return $self;
}



#Usage of new fast in the adaptor means we can't deprecate and re-assign old args in new!

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the calling code 
               is very disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::SetFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

#Remove and use Bio::EnsEMBL::Feature::new_fast? - This 'weakens' adaptor

sub new_fast {
  return bless ($_[1], $_[0]);
}



=head2 feature_set

  Example    : my $set = $efeature->feature_set();
  Description: WARNING: Can now also return ResultSet aswell as FeatureSet attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet or ResultSet
  Exceptions : None
  Caller     : General
  Status     : At Risk - marked as to be removed in v67

=cut

sub feature_set {
  return $_[0]->{'set'};
}


=head2 set

  Example    : my $fset = $set_feature->set();
  Description: Getter for the set attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet or ResultSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub set {
  return $_[0]->{'set'};
}

=head2 cell_type

  Example    : my $cell_name = $efeature->cell_type()->name();
  Description: Getter for the cell_type attribute for this feature.
               May not always be set for ExternalFeatures.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub cell_type{
	return $_[0]->set->cell_type();
}

=head2 feature_type

  Example    : my $ft_name = $efeature->feature_type()->name();
  Description: Getter for the feature_type attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub feature_type{
  my $self = shift; 
  return (defined $self->{'feature_type'}) ?  $self->{'feature_type'} : $self->set->feature_type();
}



=head2 analysis

  Example    : my $analysis = $setfeature->analysis);
  Description: Getter for the Analysis attribute for this feature.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub analysis{
  my $self = shift;
  return (defined $self->{'analysis'}) ? $self->{'analysis'} : $self->set->analysis();
}


=head2 display_label

  Example    : my $label = $feature->display_label();
  Description: Getter for the display label of this feature.
               This will most likely be over-ridden by child class
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_label{
  return $_[0]->{'display_label'};
}


1;
