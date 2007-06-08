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


=head2 fetch_all_by_FeatureSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : my @dsets = $fs_adaptopr->fetch_all_by_FeatureSet($fset);
  Description: Retrieves DataSet objects from the database based on the FeatureSet.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if arg is not a valid FeatureSet
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureSet {
    my $self = shift;
    my $fset = shift;

    
    if(! ($fset && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet") && $fset->dbID())){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::FeatureSet object");
    }
	
    my $sql = "ds.feature_set_id = '".$fset->dbID()."'";


    return $self->generic_fetch($sql);	
}


=head2 fetch_all_by_ResultSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::ResultSet
  Example    : my @dsets = $fs_adaptopr->fetch_all_by_ResultSet($rset);
  Description: Retrieves DataSet objects from the database based on the ResultSet.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if arg is not a valid ResultSet
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_ResultSet {
    my $self = shift;
    my $rset = shift;

    if(! ($rset && $rset->isa("Bio::EnsEMBL::Funcgen::ResultSet") && $rset->dbID())){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::ResultSet object");
    }
	

	#self join here to make sure we get all linked result_sets
    my $sql = 'ds.data_set_id IN (SELECT ds.data_set_id from data_set ds where result_set_id='.$rset->dbID().')';


    return $self->generic_fetch($sql);	
}



=head2 fetch_all_by_feature_type_class

  Arg [1]    : string - class of associated FeatureSet FeatureType
  Example    : my @dsets = @{$ds_adaptopr->fetch_all_by_feature_type_class('HISTONE')};
  Description: Retrieves DataSet objects from the database based on the FeatureSet FeatureType class.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if no class arg defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_feature_type_class {
    my ($self, $class, $status) = @_;
  
	throw ('Must provide a FeatureType class to retrieve DataSets') if ! defined $class;
  
	my ($constraint, @dsets);

	if($status){
	  $constraint = $self->status_to_constraint($status) if $status;
    }

	foreach my $dset(@{$self->generic_fetch($constraint)}){
	  push @dsets, $dset if $dset->feature_set->feature_type->class() eq $class;
	}

	return \@dsets;	
}

=head2 fetch_all_displayable_by_feature_type_class

  Arg [1]    : string - class of associated FeatureSet FeatureType
  Arg [2]    : string - status name e.g. DISPLAYABLE
  Example    : my @dsets = @{$ds_adaptopr->fetch_all_by_feature_type_class('HISTONE')};
  Description: Wrapper method, retrieves all displayable DataSets with given FeatureSet FeatureType class
  Returntype : Listref of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_displayable_by_feature_type_class {
    my ($self, $class) = @_;
  
	return $self->fetch_all_by_feature_type_class($class, 'DISPLAYABLE');	
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
	    ds.feature_set_id  ds.name
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

#sub _default_where_clause {
#  my $self = shift;

  #Will sthis cause problems if ds or fs is absent?
	
#  return 'ds.data_set_id = s.table_id and s.table_name="data_set"';
  #unnecessary join, need to reimplment StatusAdaptor
  #This sitll returns duplicate records matching the number of records in status
#}

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
#	return ' ORDER BY fs.feature_type_id, fs.cell_type_id'; #group on cell_type_id
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
  
  my (@data_sets, $data_set, $dbID, $rset_id, $fset_id, $fset, $rset, $name);

  my $fset_adaptor = $self->db->get_FeatureSetAdaptor();
  my $rset_adaptor = $self->db->get_ResultSetAdaptor();
  $sth->bind_columns(\$dbID, \$rset_id, \$fset_id, \$name);

  
  while ( $sth->fetch() ) {



    if($data_set && ($data_set->dbID() == $dbID)){

      if((defined $data_set->feature_set() && ($fset_id == $data_set->feature_set->dbID())) ||
		 (($fset_id == 0) && (! defined $data_set->feature_set()))){

		my $rset = $rset_adaptor->fetch_by_dbID($rset_id);

		if($rset){
		  $data_set->add_ResultSet($rset);
		}else{
		  warn "DataSet $name is linked to a missing ResultSet with dbID $rset_id\n";
		}
		  

      }
	  else{
		throw("DataSet does not yet accomodate multiple feature_sets per DataSet");
      }
    }
	else{
      push @data_sets, $data_set if($data_set);

	  #handle absent sets, dbIDs of 0
      $fset = ($fset_id) ? $fset_adaptor->fetch_by_dbID($fset_id) : undef;
	  $rset = ($rset_id) ? $rset_adaptor->fetch_by_dbID($rset_id) : undef;

      $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(
													  -DBID        => $dbID,
													  -NAME        => $name,
													  -FEATURE_SET => $fset,
													  -RESULT_SET  => $rset,
													  -ADAPTOR     => $self,
													 );
    }
  }

  #we could do the sort on cell and types here
  #we can't do a default sort as some may have feature_set or result_sets absent

  push @data_sets, $data_set if($data_set);


  #As we're not (quite) constraining how DataSets are associated
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
  
  return \@data_sets;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::OligoFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given OligoFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : Listref of stored DataSet objects
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

  my $sth = $self->prepare("INSERT INTO data_set (feature_set_id, result_set_id, name) 
                            VALUES (?, ?, ?)");
  my $sth2 = $self->prepare("INSERT INTO data_set (dbid, feature_set_id, result_set_id, name) 
                            VALUES (?, ?, ?, ?)");

  my $db = $self->db();


 FEATURE: foreach my $dset (@dsets) {

    if( ! ref $dset || ! $dset->isa('Bio::EnsEMBL::Funcgen::DataSet') ) {
      throw('Must pass a DataSet object to store');
    }

    if ( $dset->is_stored($db) ) {
      throw('DataSet [' . $dset->dbID() . '] is already stored in the database, use update_sets method to store new Result/FeatureSets in this DataSet');
    }

   
    my $fset_id = (defined $dset->feature_set()) ? $dset->feature_set->dbID() : 0;
    
	my @rsets = @{$dset->get_ResultSets()};

	if(@rsets){
	  
	  foreach my $rset (@rsets){
		
		if(! ($rset->isa("Bio::EnsEMBL::Funcgen::ResultSet") && $rset->is_stored($db))){
		  throw("All ResultSets must be stored previously") if(! $dset->feature_set->is_stored($db));
		}
		
		if(! defined $dset->dbID()){
		  $sth->bind_param(1, $fset_id || 0,         SQL_INTEGER);
		  $sth->bind_param(2, $rset->dbID() || 0,    SQL_INTEGER);
		  $sth->bind_param(3, $dset->name(),    SQL_VARCHAR);
		  $sth->execute();
		  
		  $dset->dbID( $sth->{'mysql_insertid'} );
		  $dset->adaptor($self);
		  
		}
		else{
		  $sth2->bind_param(1, $dset->dbID(),   SQL_INTEGER);
		  $sth2->bind_param(2, $fset_id || 0,        SQL_INTEGER);
		  $sth2->bind_param(3, $rset->dbID() || 0,   SQL_INTEGER);  
		  $sth2->bind_param(4, $dset->name(),   SQL_VARCHAR); 
		  $sth2->execute();
		}
	  }
	}else{#got feature_set only data set

	  if(! defined $dset->dbID()){
		$sth->bind_param(1, $fset_id,         SQL_INTEGER);
		$sth->bind_param(2, 0,                SQL_INTEGER);
		$sth->bind_param(3, $dset->name(),    SQL_VARCHAR);
		$sth->execute();
		
		$dset->dbID( $sth->{'mysql_insertid'} );
		$dset->adaptor($self);
		
	  }
	  else{
		$sth2->bind_param(1, $dset->dbID(),   SQL_INTEGER);
		$sth2->bind_param(2, $fset_id,        SQL_INTEGER);
		$sth2->bind_param(3, 0,               SQL_INTEGER);  
		$sth2->bind_param(4, $dset->name(),   SQL_VARCHAR); 
		$sth2->execute();
	  }
	}
  }
      
  return \@dsets
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


1;

