#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ResultSetAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::ResultSetAdaptor - A database adaptor for fetching and
storing ResultSet objects.  

=head1 SYNOPSIS

my $rset_adaptor = $db->get_ResultSetAdaptor();

my @rsets = @{$rset_adaptor->fetch_all_ResultSets_by_Experiment()};
my @displayable_rsets = @{$rset_adaptor->fetch_all_displayable_ResultSets()};

#Other methods?
#by FeatureType, CellType all with displayable flag?


=head1 DESCRIPTION

The ResultSetAdaptor is a database adaptor for storing and retrieving
ResultSet objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::ResultSet;
use Bio::EnsEMBL::Funcgen::ResultFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(mean median);

use vars qw(@ISA);


@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

#Generates ResultSet contains info about ResultSet content
#and actual results for channel or for chips in contig set?
#omit channel handling for now as we prolly won't ever display them
#but we might use it for running analyses and recording in result_set...change to result_group or result_analyses
#data_set!!  Then we can keep other tables names and retain ResultFeature
#and change result_feature to result_set, this makes focus of result set more accurate and ResultFeatures are lightweight result objects.
#do we need to accomodate different classes of data or multiple feature types in one set?  i.e. A combi experiment (Promot + Histone mod)?
#schema can handle this...API? ignore for now but be mindful. 
#This is subtley different to handling different experiments with different features in the same ResultSet.  
#Combi will have same sample.


#This needs one call to return all displayable sets, grouped by cell_line and ordered by FeatureType
#needs to be restricted to cell line, feature type, but these fields have to be disparate from result_feature 
#as this is only a simple linker table, and connections may not always be present
#so cell tpye and feature type constraints have to be performed on load, then can assume that associated features and results
# have same cell type/feature
#so we need to group by cell_type in sql and then order by feature_type_id in sql or rearrange in code?
#This will not know about chip sets, just that a feature set is linked to various result sets
#There fore we need to use the chip_set_id or link back to the experimental_chip chip_set_ids
#this would require a self join on experimental_chip




#Result_set_id is analagous to the chip_set key, altho' we may have NR instances of the same chip set with different analysis
#if we didn't know the sets previosuly, then we would have to alter the result_set_id retrospectively i.e. change the result_set_id.#All chips in exp to be in same set until we know sets, or all in separate set?
#Do not populate data_set until we know sets as this would cause hacky updating in data_set too.


#how are we going to accomodate a combi exp?  Promot + Histone mods?
#These would lose their exp set association, i.e. same exp & sample different exp method
#we're getting close to defining the regulon here, combined results features from the same exp
#presently want them displayed as a group but ordered appropriately
#was previously treating each feature as a separate result set


#for storing/making link we don't need the Slice context
#store should check all 
#so do we move the slice context to the object methods or make optional
#then object method can check for slice and throw or take a Slice as an optional argument
#this will enable generic set to be created to allow loading and linking of features to results
#we still need to know which feature arose from which chip!!!!  Not easy to do and may span two.
#Need to genericise this to the chip_set(or use result_set_id non unique)
#We need to disentangle setting the feature to chip/set problem from the displayable problem.
#change the way StatusAdaptor works to accomodate result_set_id:table_name:table_id, as this will define unique results
#

#can we extend this to creating skeleton result sets and loading raw results too?
#

#Result.pm should be lightweight by default to enable fast web display, do we need oligo_probe_id?


#how are we going to overcome unlinked but displayable sets?
#incomplete result_feature records will be hack to update/alter?
#could have attach_result to feature method?
#force association when loading features

=head2 fetch_all_by_Experiment_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Arg [2]    : Bio::EnsEMBL::Analysis
  Arg [3]    : optional string - chip_set name
  Example    : my @rsets = @{$rset_adaptor->fetch_all_by_Experiment_Analysis($exp, $anal)};
  Description: Retrieves a list of Bio::EnsEMBL::Funcgen::ResultSets with the given Analysis from the Experiment
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ResultSet objects
  Exceptions : Throws if Analysis or Experiment are not calid and stored
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Experiment_Analysis{
  my ($self, $exp, $analysis) = @_;

  if( !($exp && $exp->isa("Bio::EnsEMBL::Funcgen::Experiment") && $exp->dbID())){
    throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::Experiment");
  }

  if( !($analysis && $analysis->isa("Bio::EnsEMBL::Analysis") && $analysis->dbID())){
    throw("Need to pass a valid stored Bio::EnsEMBL::Analysis");
  }
  

  my $constraint = "ec.experiment_id=".$exp->dbID()." AND rs.analysis_id=".$analysis->dbID();
	
  return $self->generic_fetch($constraint);
}




=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : string - type of array (e.g. AFFY or OLIGO)
  Arg [3]    : (optional) string - logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_type($slice, 'OLIGO');
  Description: Retrieves a list of features on a given slice that are created
               by probes from the specified type of array.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if no array type is provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType {
  my ($self, $ftype) = @_;

  if( !($ftype && $ftype->isa("Bio::EnsEMBL::Funcgen::FeatureType") && $ftype->dbID())){
    throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::FeatureType");
  }
	
  my $constraint = "ec.feature_type_id =".$ftype->dbID();
	
  return $self->generic_fetch($constraint);
}
 

=head2 fetch_by_name_Analysis

  Arg [1]    : string - ResultSet name
  Arg [2]    : Bio::EnsEMBL::Funcgen::Analysis
  Example    : my $rset = $rseta->fetch_by_name($exp->name().'_IMPORT');
  Description: Retrieves a ResultSet based on the name attribute
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : Throws if no name provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name_Analysis {
  my ($self, $name, $anal) = @_;

  if( ! defined $name){
    throw('Need to pass a ResultSet name');
  }

  if(!($anal && $anal->isa('Bio::EnsEMBL::Analysis') && $anal->dbID())){
	throw('You must provide a valid, stored Bio::EnsEMBL::Analysis');
  }

	
  my $constraint = "rs.name ='${name}' AND rs.analysis_id=".$anal->dbID();
	
  return $self->generic_fetch($constraint)->[0];
}


=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _tables {
  my $self = shift;
	
  return (
	  [ 'result_set',        'rs' ],
	  [ 'chip_channel',      'cc' ],
	  [ 'experimental_chip', 'ec' ],
	  [ 'channel',           'c'  ],#This causes the N(no channelrecords) records to be returned when there is no linkable channel.

	  #we can have channel here, but only if we make the link in the default where, otherwise we'll get spurious results
	  #must also make all the fetch methods use an OR constraint dependent on the table name
	  #this would also be in default where
	 );
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
	my $self = shift;

	return qw(
		  rs.result_set_id    rs.analysis_id
		  cc.table_name       cc.chip_channel_id  
          cc.table_id         ec.feature_type_id  
		  ec.cell_type_id     rs.name
		 );

	
}

=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _default_where_clause {
  my $self = shift;
	
  #return 'rs.result_set_id = cc.result_set_id AND ((cc.table_name="experimental_chip" AND cc.table_id = ec.experimental_chip_id) OR (cc.table_name="channel" AND cc.table_id = (SELECT channel_id from channel) AND ec.experimental_chip_id=(SELECT experimental_chip_id from channel)))';

  return 'rs.result_set_id = cc.result_set_id AND ((cc.table_name="experimental_chip" AND cc.table_id = ec.experimental_chip_id) OR (cc.table_name="channel" AND cc.table_id=c.channel_id AND ec.experimental_chip_id=c.experimental_chip_id))';
  

}

=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an ORDER BY clause. Sorting by oligo_feature_id would be
			   enough to eliminate duplicates, but sorting by location might
			   make fetching features on a slice faster.
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : At Risk

=cut


#do we need this?

sub _final_clause {
  return ' GROUP by cc.chip_channel_id ORDER BY rs.result_set_id, ec.cell_type_id, ec.feature_type_id';
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;
  
  my (@rsets, $last_id, $rset, $dbid, $anal_id, $anal, $ftype, $ctype, $table_id, $name);
  my ($sql, $table_name, $cc_id, $ftype_id, $ctype_id);
  my $a_adaptor = $self->db->get_AnalysisAdaptor();
  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
  my $ct_adaptor = $self->db->get_CellTypeAdaptor();
 
  $sth->bind_columns(\$dbid, \$anal_id, \$table_name, \$cc_id, \$table_id, \$ftype_id, \$ctype_id, \$name);
  
  while ( $sth->fetch() ) {

    if(! $rset || ($rset->dbID() != $dbid)){
      
      push @rsets, $rset if $rset;
      $anal = (defined $anal_id) ? $a_adaptor->fetch_by_dbID($anal_id) : undef;
      $ftype = (defined $ftype_id) ? $ft_adaptor->fetch_by_dbID($ftype_id) : undef;
      $ctype = (defined $ctype_id) ? $ct_adaptor->fetch_by_dbID($ctype_id) : undef;
            
      $rset = Bio::EnsEMBL::Funcgen::ResultSet->new(
													-DBID         => $dbid,
													-NAME         => $name,
													-ANALYSIS     => $anal,
													-TABLE_NAME   => $table_name,
													-FEATURE_TYPE => $ftype,
													-CELL_TYPE    => $ctype,
													-ADAPTOR      => $self,
												   );
    }
    
    #This assumes logical association between chip from the same exp, confer in store method?????????????????

    if(defined $rset->feature_type()){    
      throw("ResultSet does not accomodate multiple FeatureTypes") if ($ftype_id != $rset->feature_type->dbID());
    }
    
    if(defined $rset->cell_type()){
      throw("ResultSet does not accomodate multiple CellTypes") if ($ctype_id != $rset->cell_type->dbID());
    }

    #we're not controlling ctype and ftype during creating new ResultSets to store.
    #we should change add_table_id to add_ExperimentalChip and check in that method
    
    #add just the ids here, as we're aiming at quick web display.
    $rset->add_table_id($table_id, $cc_id);
  
  }

  push @rsets, $rset if $rset;
  
  return \@rsets;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ResultSet objects
  Example    : $rsa->store(@rsets);
  Description: Stores or updates previously stored ResultSet objects in the database. 
  Returntype : None
  Exceptions : Throws if a List of ResultSet objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut

sub store{
  my ($self, @rsets) = @_;

  throw("Must provide a list of ResultSet objects") if(scalar(@rsets == 0));
  
  my (%analysis_hash);
  
  my $sth = $self->prepare('INSERT INTO result_set (analysis_id, name) VALUES (?, ?)');
  
  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();
  
 FEATURE: foreach my $rset (@rsets) {
    
    if( ! ref $rset || ! $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet') ) {
      throw('Must be an ResultSet object to store');
    }
    
        
    if ( $rset->is_stored($db) ) {
      throw('ResultSet [' . $rset->dbID() . '] is already stored in the database\nResultSetAdaptor does not yet accomodate updating ResultSets');
      #would need to retrive stored result set and update table_ids
    }

    #above does not check if it has been generated from scratch but is identical i.e. recovery.
    #Need to check table_id and analysis and that it has the correct status


    
    if ( ! defined $rset->analysis() ) {
      throw('An analysis must be attached to the ResultSet objects to be stored.');
    }
    
    # Store the analysis if it has not been stored yet
    if ( ! $rset->analysis->is_stored($db) ) {
      warn("Will this not keep storing the same analysis if we keep passing the same unstored analysis?");
      $analysis_adaptor->store( $rset->analysis() );
    }
   

    $sth->bind_param(1, $rset->analysis->dbID(), SQL_INTEGER);
    $sth->bind_param(2, $rset->name(),           SQL_VARCHAR);
    
    $sth->execute();
    
    $rset->dbID( $sth->{'mysql_insertid'} );
    $rset->adaptor($self);
    
    $self->store_chip_channels($rset);
    
  }
  
  return \@rsets;
}


=head2 store_chip_channels

  Args       : Bio::EnsEMBL::Funcgen::ResultSet
  Example    : $rsa->store_chip_channel(@rset);
  Description: Convinience methods extracted from store to allow updating of chip_channel entries 
               during inline result processing which would otherwise be troublesome due to the need
               for a chip_channel_id in the result table before the ResultSet would normally be stored
               i.e. after it has been fully populated with data.
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : Throws if a stored ResultSet object is not provided
  Caller     : General
  Status     : At Risk

=cut


sub store_chip_channels{
  my ($self, $rset) = @_;
  
  if(! ($rset && $rset->isa("Bio::EnsEMBL::Funcgen::ResultSet"))){
    throw("You must pasas a valid Bio::EnsEMBL::Funcgen::ResultSet");
  }
  
  if ( ! $rset->is_stored($self->db()) ) {
    throw('ResultSet must be stored in the database before storing chip_channel entries');
  }
  
  my $sth = $self->prepare("
		INSERT INTO chip_channel (
			result_set_id, table_id, table_name
		) VALUES (?, ?, ?)
	");
  

  #Store and set all previously unstored table_ids
  foreach my $table_id(@{$rset->table_ids()}){
    
    if(! defined $rset->get_chip_channel_id($table_id)){
      $sth->bind_param(1, $rset->dbID(),       SQL_INTEGER);
      $sth->bind_param(2, $table_id,           SQL_INTEGER);
      $sth->bind_param(3, $rset->table_name(), SQL_VARCHAR);
      
      $sth->execute();
      $rset->add_table_id($table_id,  $sth->{'mysql_insertid'});
    }
  }
  return $rset;
}

=head2 list_dbIDs

  Args       : None
  Example    : my @rsets_ids = @{$rsa->list_dbIDs()};
  Description: Gets an array of internal IDs for all ProbeFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : general
  Status     : stable

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs('result_set');
}


=head2 fetch_ResultFeatures_by_Slice_ResultSet

  Arg[0]     : Bio::EnsEMBL::Slice - Slice to retrieve results from
  Arg[1]     : Bio::EnsEMBL::Funcgen::ResultSet - ResultSet to retrieve results from
  Arg[2]     : optional string - STATUS e.g. 'DIPLAYABLE'
  Example    : my @rfeatures = @{$rsa->fetch_ResultFeatures_by_Slice_ResultSet($slice, $rset, 'DISPLAYABLE')};
  Description: Gets a list of lightweight ResultFeatures from the ResultSet and Slice passed.
               Replicates are combined using a median of biological replicates based on 
               their mean techinical replicate scores
  Returntype : List of Bio::EnsEMBL::Funcgen::ResultFeature
  Exceptions : throws if not experimental_chip ResultSet
  Caller     : general
  Status     : At risk

=cut


#Could we also have an optional net size for grouping ResultFeature into  Nbp pseudo ResultFeatures?

sub fetch_ResultFeatures_by_Slice_ResultSet{
  my ($self, $slice, $rset, $ec_status) = @_;
  
  my (@rfeatures, %biol_reps, %rep_scores, @filtered_ids);
  my ($biol_rep, $score, $start, $end, $cc_id, $old_start, $old_end);

  
  my @ids = @{$rset->table_ids()};
  #should we do some more optimisation of method here if we know about presence or lack or replicates?
  
  if($ec_status){
    @filtered_ids = @{$self->status_filter($ec_status, 'experimental_chip', @ids)};

    if(! @filtered_ids){

      warn("No ExperimentalChips have the $ec_status status, No ResultFeatures retrieved");
      return \@rfeatures;
    }
  }
  

  warn "We need to build techrep and biol rep strucs here based on cc_id";

  #basically we need to build a hash of cc_id to biolrep value
  #Then we use the biolrep as a key, and push all techrep values.
  #this can then be resolved in the method below, using biolrep rather than cc_id

  my $sql = "SELECT ec.biological_replicate, cc.chip_channel_id from experimental_chip ec, chip_channel cc 
             WHERE cc.table_name='experimental_chip'
             AND ec.experimental_id=cc.table_id
             AND cc.table_id IN(".join(', ', (map $_->dbID(), @{$rset->get_ExperimentalChips()})).")";

  
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$biol_rep, \$cc_id);
  #could maybe do a selecthashref here?

  while($sth->fetch()){
	$biol_reps{$cc_id} = $biol_rep;
  }

  
  #we don't need to account for strnadedness here as we're dealing with a double stranded feature
  #need to be mindful if we ever consider expression
  
  #my $sql = "SELECT r.score, pf.seq_region_start, pf.seq_region_end FROM result r, probe_feature pf, chip_channel cc
  #           WHERE cc.result_set_id = ".$rset->dbID()."
  #           AND cc.table_id IN (".join(' ,', @ids).")
  #           AND cc.chip_channel_id = r.chip_channel_id
  #	         AND r.probe_id=pf.probe_id
  #             AND pf.seq_region_id='".$slice->get_seq_region_id()."'
  #             AND pf.seq_region_end>='".$slice->start()."'
  #             AND pf.seq_region_start<='".$slice->end()."'
  #             ORDER by pf.seq_region_start";
  

  #we don't need X Y here, as X Y for probe will be unique for cc_id.
  #any result with the same cc_id will automatically be treated as a tech rep

  $sql = 'SELECT r.score, pf.seq_region_start, pf.seq_region_end, cc.chip_channel_id FROM result r, '.
	'probe_feature pf, chip_channel cc WHERE cc.result_set_id = '.$rset->dbID();

  $sql .= ' AND cc.table_id IN ('.join(' ,', @filtered_ids).')' if ((@filtered_ids != @ids) && $ec_status);

  $sql .= ' AND cc.chip_channel_id = r.chip_channel_id'.
	      ' AND r.probe_id=pf.probe_id'.
		  ' AND pf.seq_region_id='.$slice->get_seq_region_id().
          ' AND pf.seq_region_end>='.$slice->start().
          ' AND pf.seq_region_start<='.$slice->end().
          ' ORDER by pf.seq_region_start'; #do we need to add probe_id here as we may have probes which start at the same place

  $sth = $self->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$score, \$start, \$end, \$cc_id);
  my $position_mod = $slice->start() + 1;
  


  while ( $sth->fetch() ) {
    #we need to get best result here if start and end the same
    
    #set start end for first result
    $old_start ||= $start;
    $old_end   ||= $end; 
    
    if(($start == $old_start) && ($end == $old_end)){#First result and duplicate result for same feature
	  push @{$rep_scores{$biol_reps{$cc_id}}}, $score;
    }else{#Found new location
   
	  #store previous feature with best result from @scores
		#Do not change arg order, this is an array object!!
      push @rfeatures, Bio::EnsEMBL::Funcgen::ResultFeature->new_fast
		([($old_start - $position_mod), 
		  ($old_end - $position_mod),
		  $self->resolve_replicates_by_ResultSet(\%rep_scores, $rset)]);

	  undef %rep_scores;
	  $old_start = $start;
      $old_end = $end;
	  
      #record new score
	  @{$rep_scores{$biol_reps{$cc_id}}} = ($score);
    }
  }
  
  #store last feature  
  #Do not change arg order, this is an array object!!
  #only if found previosu results
  if($old_start){
    push @rfeatures, Bio::EnsEMBL::Funcgen::ResultFeature->new_fast
      ([($old_start - $position_mod), 
	($old_end - $position_mod),
		$self->resolve_replicates_by_ResultSet(\%rep_scores, $rset)]);

	#(scalar(@scores) == 0) ? $scores[0] : $self->_get_best_result(\@scores)]);
  }
  
  return \@rfeatures;
}



=head2 resolve_replicates_by_ResultSet

  Arg[0]     : HASHREF - chip_channel_id => @scores pairs
  #Arg[1]     : Bio::EnsEMBL::Funcgen::ResultSet - ResultSet to retrieve results from
  Example    : my @rfeatures = @{$rsa->fetch_ResultFeatures_by_Slice_ResultSet($slice, $rset, 'DISPLAYABLE')};
  Description: Gets a list of lightweight ResultFeatures from the ResultSet and Slice passed.
               Replicates are combined using a median of biological replicates based on 
               their mean techinical replicate scores
  Returntype : List of Bio::EnsEMBL::Funcgen::ResultFeature
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut


#this may be done better inline rather than sub'd as we're going to have to rebuild the duplicate data each time?
#Rset should return a hash of cc_ids keys with replicate name values.
#these can then beused to build replicate name keys, with an array technical rep score values.
#mean each value then give median of all.


sub resolve_replicates_by_ResultSet{
  my ($self, $rep_ref) = @_;#, $rset) = @_;

  my ($score, @scores, $biol_rep);

  #deal with simplest case first and fastest?
  #can we front load this with the replicate set info, i.e. if we know we only have one then we don't' have to do all this testing and can do a mean

  if(scalar(keys %{$rep_ref}) == 1){
	
	($biol_rep) = keys %{$rep_ref};

	@scores = @{$rep_ref->{$biol_rep}};

	if (scalar(@scores) == 1){
	  $score = $scores[0];
	}else{
	  $score = mean(\@scores)
	}
  }else{#deal with biol replicates

	foreach $biol_rep(keys %{$rep_ref}){
	  push @scores, mean($rep_ref->{$biol_rep});
	}

	$score = median(\@scores);

  }

  return $score;
}


sub _get_best_result{
  my ($self, $scores) = @_;


  my ($median);
  my $count = scalar(@$scores);
  my $index = $count-1;
  #need to deal with lines with no results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #deal with one score fastest
  return  $scores->[0] if ($count == 1);
  
  #taken from Statistics::Descriptive
  #remeber we're dealing with size starting with 1 but indices starting at 0

  if ($count % 2) { #odd number of scores
    $median = $scores->[($index+1)/2];
  }
  else { #even
    $median = ($scores->[($index)/2] + $scores->[($index/2)+1] ) / 2;
  }

  return $median;
}



=head2 fetch_results_by_probe_id_ResultSet

  Arg [1]    : int - probe dbID
  Arg [2]    : Bio::EnsEMBL::Funcgen::ResultSet
  Example    : my @probe_results = @{$ofa->fetch_results_by_ProbeFeature_ResultSet($pid, $result_set)};
  Description: Gets result for a given probe in a ResultSet
  Returntype : ARRAYREF
  Exceptions : throws if args not valid
  Caller     : General
  Status     : At Risk - Change to take Probe?

=cut

sub fetch_results_by_probe_id_ResultSet{
  my ($self, $probe_id, $rset) = @_;
  
  throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::ResultSet") if (! ($rset  &&
									       $rset->isa("Bio::EnsEMBL::Funcgen::ResultSet")
									       && $rset->dbID()));
  
  throw('Must pass a probe dbID') if ! $probe_id;
  
  
  
  my $cc_ids = join(',', @{$rset->chip_channel_ids()});

  my $query = "SELECT r.score from result r where r.probe_id ='${probe_id}'".
    " AND r.chip_channel_id IN (${cc_ids}) order by r.score;";

 
  #  warn "query is $query";

  my @results = map $_ = "@$_", @{$self->dbc->db_handle->selectall_arrayref($query)};
  

  return \@results;
}





1;

