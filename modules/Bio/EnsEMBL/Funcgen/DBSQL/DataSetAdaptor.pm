#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::DataSetAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::DataSetAdaptor - A database adaptor for fetching and
storing DataSet objects.  

=head1 SYNOPSIS

my $dset_adaptor = $db->get_DataSetAdaptor();

my $dset = $dset_adaptor->fetch_by_dbID(1);
my @displayable_dsets = $dset_adaptor->fetch_all_displayable_DataSets();

=head1 DESCRIPTION

The DataSetAdaptor is a database adaptor for storing and retrieving
DataSet objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);
use strict;
use warnings;

@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);#do we need this?
#we can't pass the slice through the BaseAdaptor
#So we either don't use it, or use the slice on all the DataSet calls?



#Generates DataSet contains info about DataSet content
#do we need to accomodate different classes of data or multiple feature types in one set?  i.e. A combi experiment (Promot + Histone mod)?
#schema can handle this...API? ignore for now but be mindful. 
#This is subtley different to handling different experiments with different features in the same DataSet.  
#Combi will have same sample.
#See DataSet for definitions of set types


#This needs one call to return all displayable sets, grouped by cell_line and ordered by FeatureType
#needs to be restricted to cell line, feature type, but these fields have to be disparate from result_feature 
#as this is only a simple linker table, and connections may not always be present
#so cell tpye and feature type constraints have to be performed on load, then can assume that associated features and results
# have same cell type/feature
#so we need to group by cell_type in sql and then order by feature_type_id in sql or rearrange in code?
#This will not know about chip sets, just that a feature set is linked to various result sets
#There fore we need to use the chip_set_id or link back to the experimental_chip chip_set_ids
#this would require a self join on experimental_chip



#what other methods do we need? all with displayable option
#fetch by Experiment_Slice?
#fetch by FeatureType_Slice?
#fetch by CellType_Slice?
#set_chips (duplicates)?  We are now using result_set_id as the chip_set key, so if we didn't know the sets previosuly, then we would have to alter the result_set_id retrospectively i.e. change the result_set_id.  This would also require a check on result_feature to see if any feature_sets had been associated, and cleaning up of duplicate result_feature entries if the same feature were attached to both of  the previously separate result sets.
#this may require an on duplicate key call....delete one.

#wouldn't it be better to associate the chip set info with the ec's rather than the result set?

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

=head2 fetch_all_by_Slice_ExperimentalChips

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2...] : listref of Bio::EnsEMBL::Funcgen::ExperimentalChip objects
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_arrayname($slice, $exp);
  Description: Retrieves a list of features on a given slice that are created
               by probes from the given ExperimentalChip.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoFeature objects
  Exceptions : Throws if no array name is provided
  Caller     : 
  Status     : At Risk

=cut

#This is no longer appropriate
#should this take >1 EC? What if we can't fit a all mappings onto one chip
#Would possibly miss some from the slice


sub fetch_all_by_Slice_ExperimentalChips {
	my ($self, $slice, $exp_chips) = @_;

	my (%nr);


	foreach my $ec(@$exp_chips){
				
	  throw("Need pass listref of valid Bio::EnsEMBL::Funcgen::ExperimentalChip objects") 
	    if ! $ec->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip");
	  
	  $nr{$ec->array_chip_id()} = 1;
	}

	#get array_chip_ids from all ExperimentalChips and do a
	#where op.array_chip_id IN (".(join ", ", @ac_ids)

	#my @echips = @{$self->db->get_ExperimentalChipAdaptor->fetch_all_by_experiment_dbID($exp->dbID())};
	#map $nr{$_->array_chip_id()} = 1, @echips;
	my $constraint = " op.array_chip_id IN (".join(", ", keys %nr).") AND op.oligo_probe_id = of.oligo_probe_id ";


	
	return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}




=head2 fetch_all_by_Slice_type

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

sub fetch_all_by_Slice_type {
	my ($self, $slice, $type, $logic_name) = @_;

	throw("Not implemented yet\n");
	
	throw('Need type as parameter') if !$type;
	
	my $constraint = qq( a.type = '$type' );
	
	return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint, $logic_name);
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
	  [ 'data_set',    'ds' ], 
	  [ 'result_set',  'rs' ], 
	  [ 'feature_set', 'fs' ],
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

  #will this work? May have multiple record/result_set_id
  
  
  return qw(
	    ds.data_set_id     ds.result_set_id
	    ds.feature_set_id
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
	
  return 'ds.result_set_id = rs.result_set_id AND ds.feature_set_id = fs.feature_set_id';
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


#sub _final_clause {
#	return ' ORDER BY ds.data_set_id, fs.feature_type_id, fs.cell_type_id'; #group on cell_type_id
#}


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
  
  my (@datasets, $data_set, $dbID, $rset_id, $fset_id, $fset);

  #analysis/ctype/ftype hashed
  #my (%ft_cache, %ct_chace, %anal_cache);

  my $fset_adaptor = $self->db->get_FeatureSetAdaptor();
  my $rset_adaptor = $self->db->get_ResultSetAdaptor();
  #my $anal_adaptor = $self->db->get_AnalysisAdaptor();

  #How should we store ResultSets for access? result_sets hash keyed on analysis/analysis_id?, with arrays of results sets as elements
  
  $sth->bind_columns(\$dbID, \$rset_id, \$fset_id);
  
  while ( $sth->fetch() ) {
    #This needs to be quite clever here and create a new DataSet when we encounter a new data_set_id
    #otherwise we populate the current DataSet with more result/feature sets and set them according to there analysis id?
    #feature_type?

    #do we need to check that the feature_set.sell_type_id is the same as the experimental_chip.cell_line_id


    if(! $data_set || ($data_set->dbID() == $data_set_id)){


      #we're just dealing with the basic one feature, one cell type set here.
      #Maybe we need data groups to handle anything more complex?


      #RIGHT THEN!!!
      #We need to account for non-existent feature_sets as we may only have raw data?
      $fset = (defined $fset_id) ? $fset_adaptor->fetch_by_dbID($fset_id);

      

      $data_set = Bio::EnsEMBL::Funcgen::DataSet->new( 
													  -DBID        => $dbID,
													  -FEATURE_SET => $fset,
													  -ADAPTOR     => $self,
													  #do all the rest dynamically?
													 );
    }
      #Add more result/feature sets to the dataset
      
      #Need DataSet->contains_feature_set_id($id) method
      #don't need contains_result_set_id as key will confer uniqueness on import?
      #we're likely to encounter NR feature_set_ids
      #These should be used to key which result_sets are displayed

      #As we're not constraining how DataSets are associated
      #it's valid to have a combined exp, i.e. Promoter plus Histone features and results?
      #It can be valid to have several feature_types in the same set?
      #This is entirely dependent on sensible experimental design and how we want to view the data.
      #We could then have multiple cell types too? Getting too many dimensions to display sensibly within e!


      #Logically we need one constant to key on, cell_type, feature_type, but also allow support combined info 
      #i.e. all Histone mods for one cell_type, but also have promoter data too?
      #This is getting close to a 'regulon', as we're incorporating all sorts of supporting data
      #There should be an order of display for and within each feature class 
      #(or if we have same feature across several cell types then we order alphabetically?)
      #Other possible variables to order on:  
      #analysis_id, this would also separate the features as we can't base a feature on mutliple analyses of the same data
      


      #So we either accomodate everything, where the only contraint is that we have one constant in the set
      #Or we restrict the Set to handle just one feature_set and it's supporting result_sets


    #Start simple, let's just take the one feature/data set problem first

    if($fset_id == $data_set->feature_set->dbID()){
      $data_set->add_ResultSet($rset_adaptor->fetch_by_dbID($rset_id));
    }else{
      throw("DataSet does not yet accomodate multiple feature_sets per DataSet");
    }
  }
  return \@result;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::OligoFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given OligoFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if a list of OligoFeature objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut


#We need to control what can be stored, so we need to check cell_line_ids?
#Or is this level of control to be implicit?  Would there be a use for a multi-celled DataSet
#Yes!  So we let the DB take anything, and make the obj_from_sth method handle all


sub store{
  my ($self, @dsets) = @_;

	if (scalar(@ofs) == 0) {
		throw('Must call store with a list of OligoFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO oligo_feature (
			seq_region_id,  seq_region_start,
			seq_region_end, seq_region_strand,
            coord_system_id,
			oligo_probe_id,  analysis_id,
			mismatches, cigar_line
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
	");

	my $db = $self->db();
	my $analysis_adaptor = $db->get_AnalysisAdaptor();

	FEATURE: foreach my $of (@ofs) {

		if( !ref $of || !$of->isa('Bio::EnsEMBL::Funcgen::OligoFeature') ) {
			throw('Feature must be an OligoFeature object');
		}

		if ( $of->is_stored($db) ) {
			warning('OligoFeature [' . $of->dbID() . '] is already stored in the database');
			next FEATURE;
		}

		if ( !defined $of->analysis() ) {
			throw('An analysis must be attached to the OligoFeature objects to be stored.');
		}

		# Store the analysis if it has not been stored yet
		if ( !$of->analysis->is_stored($db) ) {
			$analysis_adaptor->store( $of->analysis() );
		}

		my $original = $of;
		my $seq_region_id;
		($of, $seq_region_id) = $self->_pre_store($of);

		$sth->bind_param(1, $seq_region_id,        SQL_INTEGER);
		$sth->bind_param(2, $of->start(),          SQL_INTEGER);
		$sth->bind_param(3, $of->end(),            SQL_INTEGER);
		$sth->bind_param(4, $of->strand(),         SQL_TINYINT);
		$sth->bind_param(5, $of->coord_system_id(),SQL_INTEGER);
		$sth->bind_param(6, $of->probe->dbID(),    SQL_INTEGER);
		$sth->bind_param(7, $of->analysis->dbID(), SQL_INTEGER);
		$sth->bind_param(8, $of->mismatchcount(),  SQL_TINYINT);
		$sth->bind_param(9, $of->cigar_line(),     SQL_VARCHAR);

		$sth->execute();

		$original->dbID( $sth->{'mysql_insertid'} );
		$original->adaptor($self);
	}

	return \@ofs
}


#store_updated_sets

=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs('data_set');
}


# method by analysis and experiment/al_chip

# All the results methods may be moved to a ResultAdaptor

=head2 fetch_results_by_channel_analysis

  Arg [1]    : int - OligoProbe dbID
  Arg [2]    : int - Channel dbID
  Arg [1]    : string - Logic name of analysis
  Example    : my @results = @{$ofa->fetch_results_by_channel_analysis($op_id, $channel_id, 'RAW_VALUE')};
  Description: Gets all analysis results for probe on given channel
  Returntype : ARRAYREF
  Exceptions : warns if analysis is not valid in Channel context
  Caller     : OligoFeature
  Status     : At Risk - rename fetch_results_by_probe_channel_analysis

=cut



1;

