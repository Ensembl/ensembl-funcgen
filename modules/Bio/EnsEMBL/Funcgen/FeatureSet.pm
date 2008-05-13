#
# Ensembl module for Bio::EnsEMBL::Funcgen::FeatureSet
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::FeatureSet - A module to represent FeatureSet.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::FeatureSet;

my $result_set = Bio::EnsEMBL::Funcgen::FeatureSet->new(

); 



=head1 DESCRIPTION

A FeatureSet object provides access to a set of feature predictions and their details, which may have been generated from a 
single or multiple Experiments with potentially differing analyses.  The FeatureSet itself will only have a single analysis 
which may be one or a combination of programs but will be represented by one analysis record.


=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::FeatureSet;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate);
use Bio::EnsEMBL::Funcgen::Set;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Set);


=head2 new

  Arg [-EXPERIMENT_ID]     : Experiment dbID
  #or
  #Arg [-EXPERIMENT]       : Bio::EnsEMBL::Funcgen::Experiment
  Arg [-SLICE]             : Bio::EnsEMBL::Slice


  Example    : my $feature = Bio::EnsEMBL::Funcgen::FeatureSet->new(
                                                                    -dbid        => $dbid,
                                                                    -analysis    => $analysis,
                                                                    -feature_type => $ftype,
                                                                    -cell_type => $ctype,
                                                                    -name => $name,
                                                                    -type => 'annotated',
                                                                    -description => 'Release 3.1',
			                                                       ); 
  Description: Constructor for FeatureSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : Throws if no experiment_id defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  my ($type, $desc)
    = rearrange(['TYPE', 'DESCRIPTION'],@_);

  throw ('Must provide a FeatureType') if(! defined $self->feature_type);

  throw("You must defined a FeatureSet type e.g. 'annotated'") if ! $type;
  $self->type($type);
  $self->description($desc) if defined $desc;

  return $self;
}

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

sub new_fast {
   my ($class, $hashref)  = @_;

   return bless ($hashref, $class);
}





#methods
#set wide display label(predicted_feature) + more wordy label for wiggle tracks?
#defined by experiment type i.e. time course would require timepoint in display label
#deal with this dynamically or have display_label in table
#Need call on type, or fetch all would

#_get_ec_ids or contigsets?
#this should now be an intrinsic part of this class/adaptor

#cell line
#feature_type
#displayable...should have one for the whole set and one for each raw and predicted?

#have analysis as arg? Or do we get all analysis sets?
#we need to be able to set analyses for FeatureSets dynamically from DB
#pick up all FeatureSets 
#displayable field in FeatureSets also?

#If we have mixed types in the same experiment then we could get promoter features and histone wiggle tracks displayed togeter
#Not v.good for display purposes?  We may want to separate the promoter and histone tracks, or we may want ll the experiment data together but of mixed types.
#We need to be able to pull back the experiment type for each set, therefore this needs setting on an ec level, not an experiment level.
#This is also v.reliant on putting contig set info in place, otherwise we may get mixed chip types in same set.

#get_raw_analysis_name
#get_predicted_feature_analysis_name
#set ResultFeatures and AnnotatedFeatures in hash keyed by analysis_name?

=head2 description

  Example    : print "Feature set description is:\t".$fset->description."\n";
  Description: Getter/Setter for the description of this FeatureSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub description {
  my $self = shift;
     	
  $self->{'description'} = shift if @_;

  return $self->{'description'};
}



=head2 type

  Example    : my $type = $fset->type('annotated');
  Description: Getter/Setter for the type of this FeatureSet.
               Valid values are 'annotated', 'regulatory' or 'external'.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk 

=cut

sub type {
  my $self = shift;
     	
  #add validation of enum here..or leave to DB
  $self->{'type'} = shift if @_;

  return $self->{'type'};
}


=head2 display_label

  Example    : print $rset->display_label();
  Description: Getter for the display_label attribute for this FeatureSet.
               This is more appropriate for teh predicted_features of the set.
               Use the individual display_labels for each raw result set.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
  my $self = shift;
  
 
  if(! $self->{'display_label'}){

	if($self->feature_type->class() eq 'Regulatory Feature'){
	  $self->{'display_label'} = 'Regulatory Features';
	}
	else{
	  $self->{'display_label'} = $self->feature_type->name()." - ".$self->cell_type->name()." Enriched Sites";
	}
  }
	
  return $self->{'display_label'};
}



=head2 get_FeatureAdaptor

  Example    : 
  Description: Retrieves and caches FeatureAdaptor of feature_set type 
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::ucfirst($self->type())FeatureAdaptor
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_FeatureAdaptor{
  my $self = shift;

  if(! exists $self->{'adaptor_refs'}){

	#can we code ref this?

	$self->{'adaptor_refs'} = {(
								annotated  => $self->adaptor->db->get_AnnotatedFeatureAdaptor,
								regulatory => $self->adaptor->db->get_RegulatoryFeatureAdaptor,
								external   => $self->adaptor->db->get_ExternalFeatureAdaptor,
							   )};

  }
  return $self->{'adaptor_refs'}->{$self->type()};

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

  return $self->get_FeatureAdaptor->fetch_all_by_Slice_FeatureSets($slice, $self);
}

=head2 get_Features_by_FeatureType

  Arg[0]     : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : my @features = @{$FeatureSet->get_Features_by_FeatureType($ftype)};
  Description: Retrieves all Features for this FeatureSet for a given FeatureType
  Returntype : ARRAYRE
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_Features_by_FeatureType{
  my ($self, $type) = @_;

  return $self->get_FeatureAdaptor->fetch_all_by_FeatureType_FeatureSets($type, $self);
}





1;
