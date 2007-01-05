#
# Ensembl module for Bio::EnsEMBL::Funcgen::DataSet
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DataSet - A module to represent DataSet object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DataSet;

my $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(
	-EXPERIMENT_ID         => $exp_id,
        -SLICE                 => $slice,
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

COMPLEX/COMPOSITE - Multiple feature_types, feature classes, cell_types etc... Where the only assumtion
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
use Bio::EnsEMBL::Utils::Exception qw( throw warn );
#use Bio::EnsEMBL::Feature;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-EXPERIMENT_ID]     : Experiment dbID
  #or
  #Arg [-EXPERIMENT]       : Bio::EnsEMBL::Funcgen::Experiment
  Arg [-SLICE]             : Bio::EnsEMBL::Slice


  Example    : my $feature = Bio::EnsEMBL::Funcgen::DataSet->new(
                                                                   -EXPERIMENT_ID => $exp_id,
                                                                   -SLICE         => $slice,
                                                                   -FEATURE_TYPE  => 'HISTONE',
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
  my ($ds_id, $fs_id, $rs_id,  $ft_id, $f_anal_id, $r_anal_id, $r_table_name, $r_table_id, $cell_id)
    = rearrange(['dbID', 'FEATURE_SET_ID', 'RESULT_SET_ID', 'FEATURE_TYPE_ID', 'FEATURE_ANALYSIS_ID', 'RESULT_ANALYSIS_ID', 'RESULT_TABLE_NAME', 'RESULT_TABLE_ID', 'CELL_TYPE_ID'], @_);

 

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


  #throw("Need to pass an Experiment dbID") if ! $exp_id;
  throw("Need to pass an FeatureType dbID") if ! $ft_id;
  warn("You are defining a global context DataSet") if ! $slice;
  throw("Need to pass a CellType dbID") if ! $cell_type_id;
  
  #do we need some control of creating new objects with dbID and adding result_groups/feature_sets and them storing/updating them
  #potential for someone to create one from new using a duplicate dbID and then linking incorrect data to a pre-existing ResultGroup
  #can we check wether caller is DataSetAdaptor if we have dbID?

  if($self->dbID() && ! $caller->isa("Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor")){
    throw("You must use the DataSetAdaptor to generate DataSets with dbID i.e. from the DB, as this module accomodates updating which may cause incorrect data if the object is not generated from the DB");
  }



  warn("Need to handle single or multiple experiment and feature_group ids");


  #$self->experiment_id($exp_id) if $exp_id;
  #make these mutually exclusive?
  #$self->experiment_ids($exp_ids) if $exp_ids;

  #feature vars are not mandatory as we may just want to display raw data
  $self->feature_set_id($fg_id) if $fg_id;
  $self->feature_analysis_id($f_anal_id) if $f_anal_id;
  $self->feature_type_id($ft_id) if $ft_id,;
  $self->cell_type_id($cell_id) if $cell_id;


  if($r_anal_id || $r_table_name || $r_table_id){
    $self->add_result_set($rs_id, $r_anal_id, $r_table_name, $r_table_id);
    #This need to warn if last 3 are not defined
    #having rs_id not defined means we are storing a fresh DataSet or updating a previously stored DataSet
    #Unique combined key should prevent duplication, or should we split key and force check in adaptor? 
  }




  $self->slice($slice) if $slice;
  
  



  return $self;
}

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


#This will fail if we do a new_fast with result_set vars

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




=head2 experiment_id

  Arg [1]    : (optional) int - Experiment dbID
  Example    : $result_set->experiment_id($exp_id);
  Description: Getter and setter for the experiment_id for this DataSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub experiment_id {
    my $self = shift;
	
    $self->{'experiment_id'} = shift if @_;
		
    return $self->{'experiment_id'};
}


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
	
    if(@_){
      $self->{'experiment_ids'} = @$_[0];
    }
		
    #will this work?
    return $self->{'experiment_ids'} || [ $self->{'experiment_id'} ];
}


#add experiment method?


=head2 add_result_set

  Arg [1]    : (optional) int - result_group dbID
  Example    : $result_set->result_group_id($rg_id);
  Description: Getter and setter for the result_group_id for this DataSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub add_result_set {
    my ($self, $result_set) = @_;
	
    if (! $result_set->isa("Bio::Ensembl::Funcgen::ResultSet"){
      throw("Need to pass a valid Bio::EnsEMBL::Funcgen::ResultSet");
    }elsif(! $result_set->dbID()){
      throw("add_result_set does not yet accomodate unstored result_sets");
    }

    #we can't key on the rs_id unless we store it first
    #this is wierd, we can have an unstored data_set with stored result_sets.
    #writing classes for ResultSet will not overcome this, but we could force a ResultSet to bestored before we add it.

    #should ResultSet/Adaptor contain all the fetch_methods, and leave DataSet as a kind of organisational class as a single point of access.
    #DataSetAdaptor to perform the ordering according to feature/celltype
    #This will still not resolve the complex data sets which can be accomodated by the DB.
    #Maybe we can keep the data sets as simple as there are and confer the association by tracking back to the experiment?
    #Would there only ever be one experiment for a complex data_set?

    
    #Can have more than one experiment for a compound feature set, would we ever want to display raw data?
    #This is actually an easier problem unless we are displaying two feature types(i.e. complex and compound)


    $self->{'result_sets'} ||= {};
    $self->{'result_sets'}->{$result_set->dbID()} = $result_set;

		
    return;
}


=head2 result_set_ids

  Description: Getter for the result_set_ids for this DataSet.
           
  Returntype : LIST
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub result_set_ids {
    my $self = shift;
	
    return [keys %{$self->{'result_sets'}}];
}




=head2 feature_type_id

  Arg [1]    : (optional) int - FeatureType dbID
  Example    : $result_set->feature_type_id($ft_id);
  Description: Getter and setter for the feature_type_id for this DataSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub feature_type_id {
    my $self = shift;
	
    $self->{'feature_type_id'} = shift if @_;
		
    return $self->{'feature_type_id'};
}


=head2 feature_type

  Example    : $display_label = $rset->feature_type()->class().':'.$rset->name();
  Description: Getter and setter for the feature_type_id for this DataSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub feature_type {
    my $self = shift;
	
    if(! $self->{'feature_type'}){

      if($self->adaptor()){
	$self->{'feature_type'} = $self->adaptor->db->get_FeatureTypeAdaptor->fetch_by_dbID($self->feature_type_id());
      }else{
	throw("You need to set and adaptor to retrieve the FeatureType");
      }
    }
		
    return $self->{'feature_type'};
}


=head2 feature_set_id

  Arg [1]    : (optional) int - feature_set dbID
  Example    : $data_set->feature_set_id($fs_id);
  Description: Getter and setter for the feature_group_id for this DataSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub feature_group_id {
  my $self = shift;
	
  $self->{'feature_group_id'} = shift if @_;
		
  return $self->{'feature_group_id'};
}



=head2 cell_type_id

  Arg [1]    : (optional) int - CellLine dbID
  Example    : $data_set->cell_type_id($cell_type_id);
  Description: Getter and setter for the cell_type_id for this DataSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub cell_type_id {
    my $self = shift;
	
    $self->{'cell_type_id'} = shift if @_;
    
    return $self->{'cell_type_id'};
}


=head2 cell_type

  Example    : 
  Description: Getter and setter for the feature_type_id for this DataSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub cell_type {
  my $self = shift;
	
  if(! $self->{'cell_type'}){

    if($self->adaptor()){
      $self->{'cell_type'} = $self->adaptor->db->get_CellTypeAdaptor->fetch_by_dbID($self->cell_type_id());
    }else{
      throw("You need to set and adaptor to retrieve the CellType");
    }
  }
		
  return $self->{'cell_type'};
}


=head2 slice

  Arg [1]    : (optional) - Bio::EnsEMBL::Slice
  Example    : my $rset_slice = $rset->slice();
  Description: Getter and setter for the Slice of this DataSet.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Throws if arg is not a Slice
  Caller     : General
  Status     : At Risk

=cut

sub slice {
    my $self = shift;
	

    if(@_ && (! $_[0]->isa("Bio::EnsEMBL::Slice"))){
      throw("Arg must be a Bio::EnsEMBL::Slice");
    }


    $self->{'slice'} = shift if @_;
	
    return $self->{'slice'};
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
  
  if(! $self->{'display_label'}){
    $self->{'display_label'} = $self->feature_type->name()." -";
    $self->{'display_label'} .= " ".$self->cell_type->display_name() if $self->cell_type()->display_name();
    $self->{'display_label'} .= " Enriched Sites";
  }
	
  return $self->{'display_label'};
}


=head2 get_all_displayable_results

  Example    : my @results = @{$DataSet->get_all_displayable_results()};
  Description: wrapper to get_all_results with displayable flag passed
  Returntype : List ref to array refs containing ($display_label, @result_features);
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_all_displayable_results{
  my $self = shift;
  return $self->get_all_results(1);
}

=head2 get_all_results

  Example    : my @results = @{$DataSet->get_all_displayable_results()};
  Description: Getter and lazy loader of results attribute for
               the DataSet
  Returntype : List ref to array refs containing ($display_label, @result_features);
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_all_results{
  my ($self, $displayable) = @_;

  #do we need to restrict to Analysis? Maintain web displayable focus for now
  #what about type of display?  Wiggle vs. heatbar?
  #can we return this info, or should this be part of the query?


  if(! $self->{'results'}){

    $self->{'results'} = [];

    if($self->adaptor()){

      # we need to get all sets with ecs from the same exp with same feature_type and cell_type id
      # separate on chip_set_id and analysis?
      

      #this should really be in the EC adaptor even tho the data is in the result_set table
      #but we ideally want to update this result_set table after an import using the DataSet adaptor?
      #load is v.different to fetch as we're only retrieving start, stop and score for speed purposes.


      # if we have more than one dispalyable analysis we really only want it in the display label.
      my @eca_sets = @{$self->adaptor->get_experimental_chip_sets_analysis($self->experiment_ids(), $displayable)};
      my %anal_names;
      #do some funky map to hash thing here and check the scalar keys == 1 else set append_anal flag

      #also need to check for more than one set with the same analysis here, as these will be duplicates.
      #append display label with ec uids?

      #There is also the potential to get two identical track names from different experiments
      #this is currently only controlled by the displayable function


      foreach my $eca_set(@eca_sets){
	#my $display_label = "H3K9ac - Human Bone Osteosarcoma Epithelial Cells (U2OS)";
	my $diplay_label = $self->feature_type->name().' - '.$self->cell_type->description().' ('.$self->cell_type->display_name().')';

	if($append_anal){
	  if(! exists $anal_names{$eca_set->[1]}){
	    $anal_names{$eca_set->[1]} = $self->adaptor->db->get_AnalysisAdaptor->fetch_by_dbID($eca_set->[1])->logic_name();
	  }

	  $display_label .= ':'. $anal_names{$eca_set->[1]};
	}

	push @{$self->{'results'}},  [ $display_label, $self->adaptor->fetch_all_results_by_Slice_analysis_experimental_chips($self->slice(), $eca_set) ];
	

    }else{
      throw("Need to set an adaptor to retrieve Results");
    }
  }

  return $self->{'results'};
}





=head2 get_results_by_channel_id

  Arg [1]    : int - channel_id (mandatory)
  Arg [2]    : string - Analysis name e.g. RawValue, VSN (optional)
  Example    : my @results = $feature->results();
  Description: Getter, setter and lazy loader of results attribute for
               OligoFeature objects.
  Returntype : List ref to arrays containing ('score', 'Analysis logic_name');
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_results_by_channel_id {
    my $self = shift;
    my $channel_id = shift;
    my $anal_name = shift;

    warn "This method not fully implemented, remove/deprecate?";

    #$self->{'results'} ||= {};
    $self->{'results_complete'} ||= 0;
	
    if(! $self->{'results'} || ($anal_name && ! exists $self->{'results'}{$anal_name})){
      #fetch all, set complete set flag
      $self->{'results_complete'} ||= 1 	if(! $anal_name);
      
      foreach my $results_ref(@{$self->adaptor->fetch_results_by_channel_analysis($self->probe->dbID(), 
										  $channel_id, $anal_name)}){
	
	$self->{'results'}{$$results_ref[1]} = $$results_ref[0];
      }
    }
    
    return $self->{'results'}
}


#The experiment/al chip specificity has already been done by the ofa->fetch_all_by_Slice_Experiment
#This may be called with no preceding Experiment specificity
#this would return results for all experiments
#do we need to set a default Experiment?


#THis should return both Chip and Channel based results
#just Chip for now
#maybe retrieve and hash all if not Analysis object passed?  Then return what?  


=head2 get_result_by_Analysis_ExperimentalChips

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : listref - Bio::EnsEMBL::Funcgen::ExperimentalChip
  Example    : my $result = $feature->get_result_by_Analysis_ExperimentalChips($anal, \@echips);
  Description: Getter of results attribute for a given Analysis and set of ExperimentalChips
  Returntype : float
  Exceptions : Throws is no Analysis or ExperimentalChips are not passed?
  Caller     : General
  Status     : High Risk

=cut


#make ExperimentalChips optional?

#or have DataSetAdaptor?  Do we need a DataSet?
#may not have ExperimentalChip, so would need to return ec dbID aswell


######This will break/return anomalous if
#ECs are passed from different experiments
#ECs are passed from different Arrays


sub get_result_by_Analysis_ExperimentalChips{
    my ($self, $anal, $exp_chips) = @_;

    throw("Need to pass listref of ExperimentalChips") if(scalar(@$exp_chips) == 0);
    throw("Need to pass a valid Bio::EnsEMBL::Analysis") if ! $anal->isa("Bio::EnsEMBL::Analysis");

    my (%query_ids, %all_ids, %ac_ids);
    my $anal_name = $anal->logic_name();
    
    foreach my $ec(@$exp_chips){
				
      throw("Need to pass a listref of Bio::EnsEMBL::Funcgen::ExperimentalChip objects") 
	if ! $ec->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip");

		#my $tmp_id = $self->adaptor->db->get_OligoArrayAdaptor->fetch_by_array_chip_dbID($ec->array_chip_id())->dbID();
      
		#$array_id ||= $tmp_id;
      
      #throw("You have passed ExperimentalChips from different if($array_id != $tmp_id)
      
      #if(exists  $ac_ids{$ec->array_chip_id()}){
#	throw("Multiple chip query only works with contiguous chips within an array, rather than duplicates");
 #     }
      
      $ac_ids{$ec->array_chip_id()} = 1;
      $all_ids{$ec->dbID()} = 1;
      $query_ids{$ec->dbID()} = 1 if(! exists $self->{'results'}{$anal_name}{$ec->dbID()});
      
    }
    
    
    my @ec_ids = keys %query_ids;
    my @all_ids = keys %all_ids;
    
    
    #warn "ec ids @ec_ids\n";
    #warn "all ids @all_ids\n";
    
    #$self->{'results'} ||= {};
    #$self->{'results_complete'} ||= 0;#do we need this now?
    
    if((scalar(@all_ids) - scalar(@ec_ids))> 1){
      throw("DATA ERROR - There is more than one result stored for the following ExperimentalChip ids: @all_ids");
    }		
    elsif(! $self->{'results'} || (($anal_name && scalar(@ec_ids) > 0) && scalar(@all_ids) == scalar(@ec_ids))){
      #fetch all, set complete set flag
      #$self->{'results_complete'} ||= 1 	if(! $anal_name);
      #would need to look up chip and channel analyses here and call relevant fetch
      #or pass the chip and then build the query as = or IN dependent on context of logic name
      #if there are multiple results, last one will overwrite others
      #could do foreach here to deal with retrieving all i.e. no logic name
      #Can supply mutliple chips, but probe ids "should" be unique(in the DB at least) amongst contiguous array_chips
      #build the cache based on logic name and table_id
      #cahce key??  should we cat the ec_ids together?

      my @result_refs = @{$self->adaptor->fetch_results_by_probe_experimental_chips_analysis($self->probe->dbID(), 
											     \@ec_ids, 
											     $anal_name)};

      #Remove lines with no result
      while(@result_refs && (! $result_refs[0]->[0])){
	shift @result_refs;
      }

      my $num_results = scalar(@result_refs);
      my ($result, $mpos);
      #throw("Fetched more than one result for this OligoFeature, Analysis and ExperimentalChips") if (scalar(@result_refs) >1);

      #No sort needed as we sort in the query

      if($num_results == 1){
	$result = $result_refs[0]->[0];
      }
      elsif($num_results == 2){#mean
	$result = ($result_refs[0]->[0] + $result_refs[1]->[0])/2;
    
      }
      elsif($num_results > 2){#median or mean of median flanks
	$mpos = $num_results/2;
    
	if($mpos =~ /\./){#true median
	  $mpos =~ s/\..*//;
	  $mpos ++;
	  $result =  $result_refs[$mpos]->[0];
	}else{
	  $result = ($result_refs[$mpos]->[0] + $result_refs[($mpos+1)]->[0])/2 ;
	}
      }
      
      $self->{'results'}{$anal_name}{":".join(":", @ec_ids).":"} = $result;
    }

	#do we return the ec ids here to, or do we trust that the user will know to only pass contiguous rather than duplicate chips

	#how are we going to retrieve the result for one of many possible ec id keys?
	#options, cat ec dbids as key, and grep them to find full key, then return result
	#this may hide the duplicate chip problem
	#If a query has already been made and cached,another query with one differing ID(duplicate result) may never be queried as we already have a cahced result
	#We shoulld pick up duplicates before this happens
	#If we try and mix ExperimentalChips from different experiments, then this would also cause multiple results, and hence hide some data
	
	my @keys;
	foreach my $id(@all_ids){
	  my @tmp = grep(/:${id}:/, keys %{$self->{'results'}{$anal_name}});
	  #Hacky needs sorting, quick fix for release!!

	  if(@tmp){
	    push @keys, grep(/:${id}:/, keys %{$self->{'results'}{$anal_name}});

	    last;
	  }

	}

    throw("Got more than one key for the results cache") if scalar(@keys) > 1;

    return $self->{'results'}{$anal_name}{$keys[0]};
}


#Will this be too slow, can we not do one query across all tables


1;

