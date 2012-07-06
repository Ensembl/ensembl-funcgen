#
# Ensembl module for Bio::EnsEMBL::Funcgen::DataSet
#
# You may distribute this module under the same terms as Perl itself


=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Funcgen::DataSet - A module to represent DataSet object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DataSet;

my $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(
	                                              -DBID            => $dbID,
							 					  -ADAPTOR         => $self,
                                                  -SUPPORTING_SETS => [$rset],
                                                  -FEATURE_SET     => $fset,
                                                  -DISPLAYABLE     => 1,
                                                  -NAME            => 'DATASET1',
                                                  );



=head1 DESCRIPTION

A DataSet object provides access to either or both raw results and AnnotatedFeatures
for a given experiment within a Slice, associated with set wide experimental meta data.
This was aimed primarily at easing access to data via the web API by creating
a wrapper class with convenience methods.  The focus of this class is to contain raw and
associated processed/analysed data to be displayed as a set within the browser i.e. an 
experiment may have different cell lines, features or time points, these would require different DataSets.
# However a DataSet may contain mixed data types i.e. promoter & histone???? No give separate sets?
May have duplicates for raw data but only one predicted features track??
The data in this class is kept as lightweight as possible with data being loaded dynamically.


=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DataSet;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate);
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);
#Should not be a Set as is sufficiently different


=head2 new



  Example    : my $dset = Bio::EnsEMBL::Funcgen::DataSet->new(
                                                             -SUPPORTING_SETS => [$fset1, $fset2],
                                                             -FEATURE_SET     => $fset,
                                                             -DISPLAYABLE     => 1,
                                                             -NAME            => 'DATASET1',
			                                                 );

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
  my ($fset, $sets, $name)
    = rearrange(['FEATURE_SET', 'SUPPORTING_SETS', 'NAME'], @_);
  
  
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
    throw('You must use the DataSetAdaptor to generate DataSets with dbID i.e. from the DB,'.
          ' as this module accomodates updating which may cause incorrect data if the object'.
          ' is not generated from the DB');
  }
  
  
  $self->{'supporting_sets'} ||= {};
  #throw("Must specify at least one Result/FeatureSet") if((! $sets) && (! $fset));
  #removed this to allow generation of DataSets without feature sets
  #could reimplement this if we change the DataSetAdaptor::_obj_from_sth

  $self->add_supporting_sets($sets) if $sets;
  $self->product_FeatureSet($fset)   if $fset;	
  $self->name($name)   if $name;	
  
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
#set ResultFeatures and AnnotatedFeatures in hash keyed by analysis_name?

#Need to change to simple accessor
#or should we maintain to provide explicit method for delineating between parent and supporting FeatureSets?
#yes, and sub the feature_type/cell_type checks


=head2 product_FeatureSet

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : $data_set->product_FeatureSet($fset);
  Description: Getter and setter for the main feature_set attribute for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : Throws not a valid FeatureSet or if main feature_set has already been set.
  Caller     : General
  Status     : At Risk - change to get_product_FeatureSet

=cut

sub product_FeatureSet {
  my ($self, $fset) = @_;
  
  if($fset){
	
	if (! ($fset && ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"))){
	  throw("Need to pass a valid Bio::EnsEMBL::Funcgen::FeatureSet")
	}
	
    if(defined $self->{'feature_set'}){
      throw("The main feature_set has already been set for this DataSet, maybe you want add_SupportingSets?");
    }
	else{
	  $self->_validate_and_set_types($fset);
	  $self->{'feature_set'} = $fset;
	}
  }
	
  return $self->{'feature_set'};
}


=head2 add_supporting_sets

  Arg [1]    : Array of Bio::EnsEMBL::Feature/ResultSet object
  Example    : $dset->add_supporting_sets($rset);
  Description: Adds Result/FeatureSets to the DataSet
  Returntype : none
  Exceptions : Throws if set not valid for supporting_set type of DataSet
               Throws if supporting_sets is not an array ref
  Caller     : General
  Status     : At Risk

=cut


sub add_supporting_sets {
  my ($self, $sets) = @_;
	
  #should we handle displayable here, and propogate to the ResultSet if update_status is set
  #is there scope to write a Funcgen::Storable, which provides convenience methods to StatusAdaptor?
  #would have to make sure Feature object also inherited from Funcgen::Storable aswell as BaseFeature

  throw("Supporting sets need to be a reference to an ARRAY:\t".$sets) if ref($sets) ne 'ARRAY';

  foreach my $set(@$sets){

	if(!(ref($set) &&  $set->isa('Bio::EnsEMBL::Funcgen::Set') && $set->set_type ne 'data' && $set->dbID)){
	  throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::Set which is not a DataSet:\t$set");
	}
	#set type cannot be data at present as it does not inherit from Set.pm



	#Only validate if we are dealing with result type data
	#As we can have various cell/feature_types for compound analyses e.g. RegulatoryFeatures

	$self->_validate_and_set_types($set) if $set->set_type() ne 'feature';
	
	#should ResultSet/Adaptor contain all the fetch_methods, and leave DataSet as a kind of organisational class as a single point of access.
	#DataSetAdaptor to perform the ordering according to feature/celltype
	#This will still not resolve the complex data sets which can be accomodated by the DB.
	#Maybe we can keep the data sets as simple as there are and confer the association by tracking back to the experiment?
	#Would there only ever be one experiment for a complex data_set?
	
	
	#Can have more than one experiment for a compound feature set, would we ever want to display raw data?
	#This is actually an easier problem unless we are displaying two feature types(i.e. complex and compound)


	$self->{'supporting_sets'}->{$set->analysis->dbID()} ||= ();
	push @{$self->{'supporting_sets'}->{$set->analysis->dbID()}}, $set;
  }
		
  return;
}


=head2 _validate_and_set_types

  Arg [1]    : Bio::EnsEMBL::Feature/ResultSet object
  Example    : $dset->_validate_and_set_types($rset);
  Description: Validates and sets DataSet cell and feature types
  Returntype : none
  Exceptions : Throws if types not valid
  Caller     : General
  Status     : At Risk

=cut


sub _validate_and_set_types{
  my ($self, $set) = @_;

  #slightly dodgy bypassing methods, but extendable

  #This currently restricts all set types to one cell and feature type
  #this is incorrect for feature_set types as we want to munge several feature and possibly cell types 
  #into one combined data set.
  #this should set it to the FeatureSet type if is feature_set data_set
  #this only works as we only validate supporting_sets if type is not feature

  for my $type('feature_type', 'cell_type'){

	if(defined $self->{$type}){
	  
	  #Need to test isa here?  Why is this passing the defined test if not set?
	  if($set->{$type}->name() ne $self->{$type}->name()){

		throw(ref($set)." $type(".$set->{$type}->name().
			  ") does not match DataSet $type(".$self->{$type}->name().")");
		
	  }
	}
	else{
	  $self->{$type} = $set->{$type};
	}
  }

  return;
}



=head2 get_supporting_sets_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen:Analysis
  Arg [2]    : (optional) status - e.g 'DISPLAYABLE'
  Example    : my $anal_sets = @{$result_set->get_ResultSets_by_Analysis($analysis)};
  Description: Getter for the SupportingSet objects of a given Analysis.
  Returntype : ARRAYREF
  Exceptions : Throws if arg is not a valid stored Bio::EnsEMBL::Anaylsis
  Caller     : General
  Status     : At Risk

=cut

sub get_supporting_sets_by_Analysis {
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
  
  foreach my $anal_rset(@{$self->{'supporting_sets'}->{$analysis->dbID()}}){
	  
	  if(! defined $status){
		  push @rsets, $anal_rset;
	  }
	  elsif($anal_rset->has_status($status)){
		  push @rsets, $anal_rset;
	  }
  }
		
  return \@rsets;
}



=head2 get_supporting_sets

  Arg [1]    : (optional) status - e.g 'DISPLAYABLE'
  Example    : my @status_sets = @{$data_set->get_supporting_sets($status)};
  Description: Getter for the ResultSets for this DataSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_supporting_sets{
  my ($self, $status, $set_type)  = @_;
  #swap the args here

  #Add analysis here and make above method wrapper

  #Validate type here
  if($set_type &&
	 ($set_type ne 'result' &&
	  $set_type ne 'feature' &&
	  $set_type ne 'input')){
	throw("You have specified an invalid supporting set type:\t$set_type");
  }


  my @ssets;

  foreach my $anal_id(keys %{$self->{'supporting_sets'}}){
    
    foreach my $sset(@{$self->{'supporting_sets'}->{$anal_id}}){

	  if(defined $status && 
		 (! $sset->has_status($status))){
		next;
	  }

	  if(defined $set_type &&
		 ($sset->set_type ne $set_type)){
		next;
	  }

	  push @ssets, $sset;
    }
  }

  return \@ssets;
}




=head2 get_displayable_supporting_sets

  Example    : my @displayable_rsets = @{$result_set->get_displayable_supporting_sets()};
  Description: Convenience method for web display
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_displayable_supporting_sets{
  my ($self, $set_type) = @_;

  return $self->get_supporting_sets('DISPLAYABLE', $set_type);
}



=head2 get_displayable_product_FeatureSet

  Example    : my $fset = $data_set->get_displayable_product_FeatureSet();
  Description: Convenience method for web display
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_displayable_product_FeatureSet{
  my $self = shift;

  return  $self->product_FeatureSet->has_status('DISPLAYABLE') ?  $self->product_FeatureSet() : undef;
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




#The following attributes are generated dynamically from the consituent Result/FeatureSets

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

	#This does not account for DataSet without a product FeatureSet
	my $fset = $self->product_FeatureSet;

	if($fset && ($fset->feature_type->class() eq 'Regulatory Feature')){
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


#sub get_type_config{
#  my ($self) = @_;
#
#  if (! defined $self->{type_config}){
#	$self->{type_config} = $self->adaptor->fetch_type_config_by_DataSet($self);
#  }
#
#  return $self->{type_config};
#}



1;

