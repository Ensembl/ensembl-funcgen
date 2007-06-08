#
# Ensembl module for Bio::EnsEMBL::Funcgen::DataSet
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DataSet - A module to represent DataSet object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DataSet;

my $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(
	                           -DBID        => $dbID,
							   -FEATURE_SET => $fset,
							   -ADAPTOR     => $self,
); 



=head1 DESCRIPTION

A DataSet object provides access to either or both raw results and PredictedFeatures
for a given experiment within a Slice, associated with set wide experimental meta data.
This was aimed primarily at easing access to data via the web API by creating
a wrapper class with convenience methods.  The focus of this class is to contain raw and
associated processed/analysed data to be displayed as a set within the browser i.e. an 
experiment may have different cell lines, features or time points, these would require different DataSets.
#However a DataSet may contain mixed data types i.e. promoter & histone???? No give separate sets?
May have duplicates for raw data but only one predicted features track??
The data in this class is kept as lightweight as possible with data being loaded dynamically.


SOME IMPORTANT ISSUES DEFINITIONS

This class current only accomodates the following relationships:


SIMPLE - feature_set to result_set(s) relationships.  This is one feature_set/type with a supporting 
result_set or sets from the same experiment.

COMPOUND - feature_set to result_sets relationship.  Where we have one feature_set/type supported by
numerous result_sets which may have different analyses from different experiments.


Both SIMPLE and COMPOUND also assume all other variables are the same e.g. cell_type, time_point etc.

This class does not accomodate the following:

COMPLEX - Multiple feature_types, feature classes, cell_types etc... Where the only assumtion
is that their is one constant variable which can be keyed on.  This could potentially capture any experiment design.

e.g. A combined promoter and histone tiling experiment which has features and results for promoter and all modifications,
but using the same cell line and conditions.  

e.g. Looking at the same histone modifications across multiple cell_types

e.g. Looking at time points within an experiment


Final goal of visualisation will be a track of regulons/functional features supported by a network of
feature_types/classes from different cell_types, some relationships may be indirect. 


=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DataSet;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new



  Example    : my $feature = Bio::EnsEMBL::Funcgen::DataSet->new(
                                                               -RESULT_SET     => $rset,
                                                               -FEATURE_SET    => $fset,
                                                               -DISPLAYABLE    => 1,
                                                               -UPDATE_SUBSETS => 1,
			                                        );

#for COMPLEX DataSet could use this, where 1 and 2 are the positions they are to be returned in
#Would also need to record what the display type would be for each set, so the webcode can do it dynamically.
#This would allow any config of display based on what is defined in the DB.


  #-DATA_GROUPS    => \([1, $fset1, @rsets1], [2, $fset2, @rsets2]),
  #add_ResultSet($rset, $display_priority)?

  Description: Constructor for DataSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : Throws if no experiment_id defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  #do we need to add $fg_ids to this?  Currently maintaining one feature_group focus.(combi exps?)
  my ($fset, $rset, $name)
    = rearrange(['FEATURE_SET', 'RESULT_SET', 'NAME'], @_); #, 'UPDATE_SUBSETS', 'DISPLAYABLE'], @_);

 
  my @caller = caller();

  #do we need to passexperiment_id to check that table_name/id correspond for storage?
  #'EXPERIMENT_ID', 'EXPERIMENT_IDS',

  #Can have more than one experiment_id for a combined feature set. But shouldn't query like that.
  #therefore we need to be able to track back from feature to ec's rather than exps.
  #as there may be mixed data in an exp which didn't necessarily contribute to the combined feature
  #We are now separating potentially different featuretype from the same exp into different result_groups
  #therefore we only have to track back to the result_group e.g. the contig chip set

  #We also need a way of pulling back GOLDEN/combined resultssets based on feature_set_id
  #Set status as GOLDEN, then pull back displayable or GOLDEN raw results

  #Could link experiment_feature_type/feature_set to ec or result_set table?
  #latter would mean we don't have to specifiy which ec, just part of set.
  #This will make it easier for populating pfs but will mean that we can't easily track back to a particular ec without doing some probe/slice look up via the array chip.
  #Not really a requirement, so let's take this hit.

  #Could then maybe use DataSet to store pfs, otherwise we'd have to pass the rset or at the very least the result_set_id.  
  #do we need some control of creating new objects with dbID and adding result_groups/feature_sets and them storing/updating them
  #potential for someone to create one from new using a duplicate dbID and then linking incorrect data to a pre-existing ResultGroup
  #can we check wether caller is DataSetAdaptor if we have dbID?

  if($self->dbID() && $caller[0] ne "Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor"){
    throw("You must use the DataSetAdaptor to generate DataSets with dbID i.e. from the DB, as this module accomodates updating which may cause incorrect data if the object is not generated from the DB");
  }

  #warn("Need to handle single or multiple experiment and feature_group ids");

  $self->{'result_sets'} ||= {};


  #change this to add_FeatureSet and add_ResultSet
  #both need to check whether feature or cell predefined.
  #then check names
  #add_featureSet must thro if already defined.

  throw("Must specify at least one Result/FeatureSet") if((! $rset) && (! $fset));

  $self->add_ResultSet($rset) if $rset;
  $self->feature_set($fset)   if $fset;	
  $self->name($name)   if $name;	
  
  

  #Now we need a store_sets method in DataSet Adaptor



  return $self;
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
#we need to be able to set analyses for DataSets dynamically from DB
#pick up all DataSets 
#displayable field in DataSets also?

#If we have mixed types in the same experiment then we could get promoter features and histone wiggle tracks displayed togeter
#Not v.good for display purposes?  We may want to separate the promoter and histone tracks, or we may want ll the experiment data together but of mixed types.
#We need to be able to pull back the experiment type for each set, therefore this needs setting on an ec level, not an experiment level.
#This is also v.reliant on putting contig set info in place, otherwise we may get mixed chip types in same set.

#get_raw_analysis_name
#get_predicted_feature_analysis_name
#set ResultFeatures and PredictedFeatures in hash keyed by analysis_name?



=head2 experiment_ids

  Arg [1]    : (optional) array ref - Experiment dbIDs
  Example    : $result_set->experiment_ids(\@exp_ids);
  Description: Getter and setter for the experiment_ids for this DataSet.
               Only used if the predicted features are a composite of different experiments
  Returntype : LIST
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub experiment_ids {
    my $self = shift;
	
    throw("Not yet implemented");
    #need to look up through ResultSets
		
    #will this work?
    return $self->{'experiment_ids'};
}

=head2 feature_set

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : $data_set->feature_set($fset);
  Description: Getter and setter for the feature_set attribute for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgne::FeatureSet
  Exceptions : Throws if does not match DataSet feature_type or if not a valid stored FeatureSet
  Caller     : General
  Status     : At Risk

=cut

sub feature_set {
  my ($self, $fset) = @_;
	
  #should we handle displayable here, and propogate to the FeatureSet if update_status is set
  #is there scope to write a Funcgen::Storable, which provides convenience methods to StatusAdaptor?
  #would have to make sure Feature object also inherited from Funcgen::Storable aswell as BaseFeature

  if($fset){
	
    if(defined $self->{'feature_set'}){
      throw("DataSet does not yet aqccomodate multiple FeatureSets");

    }else{

      throw("Need to pass a valid Bio::EnsEMBL::Funcgen::FeatureSet") if (! $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"));

      if(defined $self->{'feature_type'}){

	if($fset->feature_type->name() ne $self->feature_type->name()){
	  throw("FeatureSet feature_type does not match DataSet feature_type");
	}
      }else{
	$self->{'feature_type'} = $fset->feature_type();
      }

      if(defined $self->{'cell_type'}){
      
	if($fset->cell_type->name() ne $self->cell_type->name()){
	  throw("FeatureSet cell_type does not match DataSet cell_type");
	}
      }else{
	$self->{'feature_type'} = $fset->feature_type();
      }
    }

    $self->{'feature_set'} = $fset;
  }
		
  return $self->{'feature_set'};
}




=head2 add_ResultSet

  Arg [1]    : Bio::EnsEMBL::ResultSet
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : $dset->add_ResultSet($rset);
  Description: Adds ResultSets to the DataSet
  Returntype : none
  Exceptions : Throws if CellType or FeatureType do not match
  Caller     : General
  Status     : At Risk

=cut

sub add_ResultSet {
  my ($self, $rset, $displayable) = @_;
	
  #should we handle displayable here, and propogate to the ResultSet if update_status is set
  #is there scope to write a Funcgen::Storable, which provides convenience methods to StatusAdaptor?
  #would have to make sure Feature object also inherited from Funcgen::Storable aswell as BaseFeature

  if (! ($rset && $rset->isa("Bio::EnsEMBL::Funcgen::ResultSet"))){
    throw("Need to pass a valid Bio::EnsEMBL::Funcgen::ResultSet");
  }

  if(defined $self->{'feature_type'}){

    if($rset->feature_type()->name() ne $self->feature_type()->name()){
      throw("ResultSet feature_type(".$rset->feature_type->name().
	    ") does not match DataSet feature_type(".$self->feature_type->name().")");
    }
  }else{
    $self->{'feature_type'} = $rset->feature_type();
  }
  

  if(defined $self->{'cell_type'}){
  
    if($rset->cell_type()->name ne $self->cell_type()->name){
      throw("ResultSet cell_type does not match DataSet cell_type");
    }
  }else{
    $self->{'cell_type'} = $rset->cell_type();
  }
  
  
  #should ResultSet/Adaptor contain all the fetch_methods, and leave DataSet as a kind of organisational class as a single point of access.
  #DataSetAdaptor to perform the ordering according to feature/celltype
  #This will still not resolve the complex data sets which can be accomodated by the DB.
  #Maybe we can keep the data sets as simple as there are and confer the association by tracking back to the experiment?
  #Would there only ever be one experiment for a complex data_set?
  
  
  #Can have more than one experiment for a compound feature set, would we ever want to display raw data?
  #This is actually an easier problem unless we are displaying two feature types(i.e. complex and compound)
  
  $self->{'result_sets'}->{$rset->analysis->dbID()} ||= ();
  push @{$self->{'result_sets'}->{$rset->analysis->dbID()}}, $rset;

		
  return;
}


=head2 get_ResultSets_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen:Analysis
  Arg [2]    : (optional) status - e.g 'DISPLAYABLE'
  Example    : my $anal_sets = @{$result_set->get_ResultSets_by_Analysis($analysis)};
  Description: Getter for the ResultSet of given Analysis for this DataSet.
  Returntype : Arrayref
  Exceptions : Throws if arg is not a valid stored Bio::EnsEMBL::Anaylsis
  Caller     : General
  Status     : At Risk

=cut

sub get_ResultSets_by_Analysis {
  my ($self, $analysis, $status) = @_;

  my @rsets;
	

  #should we handle displayable here, and propogate to the ResultSet if update_status is set
  #is there scope to write a Funcgen::Storable, which provides convenience methods to StatusAdaptor?
  #would have to make sure Feature object also inherited from Funcgen::Storable aswell as BaseFeature

  
  if (! ($analysis->isa("Bio::EnsEMBL::Analysis") && $analysis->dbID())){
	  throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::ResultSet");
  }

  #will have to generate new array of object here if we want to filter displayable
  #This may result in returning a ref to the stored ResultSets for no status
  #And a ref to the abstracted/filtered i.e. non-stored ResultSets if we have a status
  #This could cause problems if people want to edit the real ResultSets via the refs
  #If we edit the ResultSets like this, we would still store via their adaptor
  #so would need to refresh DataSet anyway.  
  
  #should ResultSet/Adaptor contain all the fetch_methods, and leave DataSet as a kind of organisational class as a single point of access.
  #DataSetAdaptor to perform the ordering according to feature/celltype
  #This will still not resolve the complex data sets which can be accomodated by the DB.
  #Maybe we can keep the data sets as simple as there are and confer the association by tracking back to the experiment?
  #Would there only ever be one experiment for a complex data_set?
  
  
  #Can have more than one experiment for a compound feature set, would we ever want to display raw data?
  #This is actually an easier problem unless we are displaying two feature types(i.e. complex and compound)

  #could we have >1 rset with the same analysis?
  
  foreach my $anal_rset(@{$self->{'result_sets'}->{$analysis->dbID()}}){
	  
	  if(! defined $status){
		  push @rsets, $anal_rset;
	  }
	  elsif($anal_rset->has_status($status)){
		  push @rsets, $anal_rset;
	  }
  }
		
  return \@rsets;
}


=head2 get_ResultSets

  Arg [1]    : (optional) status - e.g 'DISPLAYABLE'
  Example    : my @status_sets = @{$result_set->get_ResultSets($status)};
  Description: Getter for the ResultSets for this DataSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_ResultSets{
  my ($self, $status)  = @_;

  my @rsets;

  foreach my $anal_id(keys %{$self->{'result_sets'}}){
    
    foreach my $rset(@{$self->{'result_sets'}->{$anal_id}}){

      if(! defined $status){
	push @rsets , $rset;
      }elsif($rset->has_status($status)){
	push @rsets, $rset;
      }
    }
  }

  return \@rsets;
}


=head2 get_displayable_ResultSets

  Example    : my @displayable_rsets = @{$result_set->get_displayable_ResultSets()};
  Description: Convenience method for web display
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_displayable_ResultSets{
  my $self = shift;

  return $self->get_ResultSets('DISPLAYABLE');
}


=head2 get_displayable_FeatureSets

  Example    : my @displayable_fsets = @{$result_set->get_displayable_FeatureSets()};
  Description: Convenience method for web display
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_displayable_FeatureSets{
  my $self = shift;

  #need to write get_feature_sets, only when we accomodate multiple feature_sets
  #this is just a place holder method to reduce change in teh AI with repsect to the web API

  my @fsets = ();
  
  push @fsets, $self->feature_set() if $self->feature_set->has_status('DISPLAYABLE');


  return \@fsets;
}


=head2 result_set_ids

  Description: Getter for the result_set_ids for this DataSet.
           
  Returntype : LIST
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

#sub result_set_ids {
#    my $self = shift;
	
    #return [keys %{$self->{'result_sets'}}];
#}




=head2 cell_type

  Example    : my $dset_ctype_name = $dset->cell_type->name();
  Description: Getter for the cell_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub cell_type {
  my $self = shift;
     		
  return $self->{'cell_type'};
}

=head2 feature_type

  Example    : my $dset_ftype_name = $dset->feature_type->name();
  Description: Getter for the feature_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub feature_type {
  my $self = shift;
     		
  return $self->{'feature_type'};
}

=head2 name

  Example    : my $dset->name('DATASET1');
  Description: Getter/Setter for the name of this DataSet.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub name {
  my $self = shift;
     	
  $self->{'name'} = shift if @_;

  return $self->{'name'};
}



=head2 display_label

  Example    : print $rset->display_label();
  Description: Getter for the display_label attribute for this DataSet.
               This is more appropriate for teh predicted_features of the set.
               Use the individual display_labels for each raw result set.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
  my $self = shift;
  

  #Add display label in table?

  if(! $self->{'display_label'}){

	if($self->feature_set->feature_type->class() eq 'REGULATORY FEATURE'){
	  $self->{'display_label'} = 'Regulatory Features';
	}
	else{

	  $self->{'display_label'} = $self->feature_type->name()." -";
	  $self->{'display_label'} .= " ".($self->cell_type->display_label() || 
									   $self->cell_type->description()   ||
									   $self->cell_type()->name());
	  $self->{'display_label'} .= " Enriched Sites";
	}
  }
 
  return $self->{'display_label'};
}





1;

