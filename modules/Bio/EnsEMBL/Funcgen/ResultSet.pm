#
# Ensembl module for Bio::EnsEMBL::Funcgen::ResultSet
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::ResultSet - A module to represent ResultSet.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ResultSet;

my $result_set = Bio::EnsEMBL::Funcgen::ResultSet->new(

); 



=head1 DESCRIPTION

A ResultSet object provides access to a set raw results from an Experiment. A set will be one or more 
contiguous chips to be treated as one set, with the same analysis. Duplicate sets will form a separate
result set, as will the same raw data analysed or normalised in a different manner.


=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::ResultSet;

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


  Example    : my $feature = Bio::EnsEMBL::Funcgen::ResultSet->new(
                                                                   -dbid        => $dbid,
                                                                   -analysis    => $analysis,
                                                                   -table_name  => 'experimental_chip',
                                                                   -table_id    => $ec_id,
			                                          ); 
  Description: Constructor for ResultSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : Throws if no experiment_id defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  my ($anal_id, $table_name, $table_id)
    = rearrange(['ANALYSIS', 'TABLE_NAME', 'TABLE_ID'], @_);

  #maybe don't need tha analysis args as mandatory as we're testing in the adaptor store method
  if (! ( $table_name && $table_id)){
    throw("Need to pass the following args:\ttable_name\ttable_id");
  }

 
  
  #do we need some control of creating new objects with dbID and adding result_groups/feature_sets and them storing/updating them
  #potential for someone to create one from new using a duplicate dbID and then linking incorrect data to a pre-existing ResultGroup
  #we need to verify that each table_name/id in the set is from the same experiment


  if($self->dbID() && ! $caller->isa("Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor")){
    warn("You may be adding ${table_name}:${table_id} to a previously existing ResultSet");
    #This is only true if the dbID passed has been used before
  }

  $self->analysis($analysis) if $analysis;
  $self->table_name($table_name);
  $self->add_table_id($table_id);


  #$self->experiment_id($exp_id) if $exp_id;#should have this method but only as a getter
  #$self->slice($slice) if $slice;#should always pass slice as arg as we'll only ever do it once?
  

  return $self;
}

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
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
#we need to be able to set analyses for ResultSets dynamically from DB
#pick up all ResultSets 
#displayable field in ResultSets also?

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
  Description: Getter and setter for the experiment_id for this ResultSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub experiment_id {
    my $self = shift;
	

    throw("Not yet implemented");

    #if(! defined $self->{'experiment_id'}){
      

    #$self->{'experiment_id'} = $self->adaptor->db->get_ExperimentalChipAdaptor->fetch_by_dbID(;
		
    return $self->{'experiment_id'};
}

=head2 analysis

  Arg [1]    : (optional) - Bio::EnsEMBL::Analysis
  Example    : $anal_id = $rset->analysis->dbID();
  Description: Getter and setter for the analysis attribute for this ResultSet.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub analysis {
  my $self = shift;
	
  $self->{'analysis'} = shift if(@_);
  #checking for Analysis done in store
		
  return $self->{'analysis'};
}



=head2 add_table_id

  Example    : $result_set->add_table_id($ec_id);
  Description: Appends a table_id to the ResultSet.
  Returntype : None
  Exceptions : Throws if no table_id defined
  Caller     : General
  Status     : At Risk

=cut

sub add_table_id {
  my ($self, $table_id) = @_;

  if (! defined $table_id){	
    throw("Need to pass a table_id");
  }else{
    $self->{'table_ids'} ||= [];
    push @{$self->{'table_ids'}}, $table_id;
  }

  return;
}


=head2 table_ids

  Example    : $result_set->feature_group_id($fg_id);
  Description: Getter and setter for the feature_group_id for this ResultSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub table_ids {
  my $self = shift;
			
  return $self->{'table_ids'};
}



=head2 display_label

  Example    : print $rset->display_label();
  Description: Getter for the display_label attribute for this ResultSet.
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


=head2 get_displayable_ResultFeatures_by_Slice

  Example    : my @results = @{$ResultSet->get_all_displayable_results()};
  Description: wrapper to get_all_results with displayable flag passed
  Returntype : List ref to array refs containing ($display_label, @result_features);
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_displayable_ResultFeatures_by_Slice{
  my ($self, $slice) = @_;
  return $self->get_ResultFeatures_by_Slice($slice, 1);
}

=head2 get_ResultFeatures_by_Slice

  Example    : my @results = @{$ResultSet->get_all_displayable_results()};
  Description: Getter and lazy loader of results attribute for
               the ResultSet
  Returntype : List ref to array refs containing ($display_label, @result_features);
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_ResultFeatures_by_Sliceblaart{
  my ($self, $slice, $displayable) = @_;

  #do we need to restrict to Analysis? Maintain web displayable focus for now
  #what about type of display?  Wiggle vs. heatbar?
  #can we return this info, or should this be part of the query?

  #Does this also need to accomodate channel level data?
  #No! Normalisation is never done in a slice context, rather a chip context.



  if(! $self->{'result_features'}){

    $self->{'result_features'} = [];

    if($self->adaptor()){

      # we need to get all sets with ecs from the same exp with same feature_type and cell_type id
      # separate on chip_set_id and analysis?
      

      #this should really be in the EC adaptor even tho the data is in the result_set table
      #but we ideally want to update this result_set table after an import using the ResultSet adaptor?
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


	#should do a _new_fast on ResultFeature here
	#No need for a ResultFeatureAdaptor as they are transient i.e. not storable and only access through a result set
	#this needs totally changing to be focused on one resultset



	#push @{$self->{'results'}},  [ $display_label, $self->adaptor->fetch_all_results_by_Slice_analysis_experimental_chips($self->slice(), $eca_set) ];
	@{$self->{'result_features'}} = $self->adaptor->fetch_all_results_by_Slice_analysis_experimental_chips($slice(), $self->analysis_id, $eca_set) ];
	

    }else{
      throw("Need to set an adaptor to retrieve Results");
    }
  }

  return $self->{'results'};
}



#Is it possible to to have one chip displayable and another not within the same set.
#one may fail validation, but would just omit from result_set?
#Or maybe we only want certain chips displayed, locational context is not split evenly over chips so not that useful :?
#If we handle this here then all we have to do is filter table_ids first before passing them to the fetch.

sub get_ResultFeatures_by_Slice{
  my ($self, $slice, $displayable) = @_;


  my (@ofs, @results, $result);

  
  foreach my $of(@{$self->fetch_all_by_Slice_ExperimentalChips($slice, $exp_chips)}){
    
    if((! @ofs) || ($of->start == $ofs[0]->start() && $of->end == $ofs[0]->end())){
      push @ofs, $of;
    }else{#Found new location, deal with previous
      push @results, [$ofs[0]->start(), $ofs[0]->end(), $self->_get_best_result(\@ofs, $analysis, $exp_chips)];
      @ofs = ($of);
    }
  }

  push @results, [$ofs[0]->start(), $ofs[0]->end(), $self->_get_best_result(\@ofs, $analysis, $exp_chips)];

  return \@results;

}

sub _get_best_result{
  my ($self, $ofs, $analysis, $exp_chips) = @_;

  my ($result, $mpos);

  if(scalar(@$ofs) == 2){#mean
    $result = ($ofs->[0]->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips) + 
	       $ofs->[1]->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips))/2;
    
  }
  elsif(scalar(@$ofs) > 2){#median or mean of median flanks
    $mpos = (scalar(@$ofs))/2;
    
    if($mpos =~ /\./){#true median
      $mpos =~ s/\..*//;
      $mpos ++;
      $result = $ofs->[$mpos]->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips);
    }else{
      $result = ($ofs->[$mpos]->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips) +
		 $ofs->[($mpos+1)]->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips))/2 ;
    }
  }else{
    #push start, end, score onto results
    $result =  $ofs->[0]->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips);

  }

  return $result;
}





1;
