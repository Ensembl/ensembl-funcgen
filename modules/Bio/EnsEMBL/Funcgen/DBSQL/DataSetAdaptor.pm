#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::DataSetAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::DataSetAdaptor - A database adaptor for fetching and
storing DataSet objects.  

=head1 SYNOPSIS

my $dset_adaptor = $db->get_DataSetAdaptor();

my $dset = $dset_adaptor->fetch_by_dbID(1);
my @displayable_dsets = $dset_adaptor->fetch_all_displayable();

=head1 DESCRIPTION

The DataSetAdaptor is a database adaptor for storing and retrieving
DataSet objects.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);
use strict;
use warnings;

@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor); #do we need this?
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

=head2 fetch_by_name

  Arg [1]    : string - name of DataSet
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my $dset = $dset_adaptor->fetch_by_name('RegulatoryFeatures:MultiCell');
  Description: Fetch DataSet with a given name
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : Throws if no name passed 
  Caller     : General
  Status     : At Risk 

=cut

sub fetch_by_name {
  my ($self, $name, $status) = @_;
  
  throw("Must provide a name argument") if (! defined $name);
  
  my $sql = "ds.name='".$name."'";
  
  if ($status) {
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return (defined $sql) ? $self->generic_fetch($sql)->[0] : [];
  
}

=head2 fetch_all_by_supporting_set_type

  Arg [1]    : string - type of supporting_sets i.e. result or feature
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my $dsets = $dset_adaptor->fetch_all_by_supporting_set('feature');
  Description: Fetch all DataSets whose pre-processed data consists of a particular set type
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if no supporting_set_type passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_supporting_set_type {
  my ($self, $type, $status) = @_;
  
  throw("Must provide a supporting_set type argument") if (! defined $type);
  
  my $sql = "ss.type='".$type."'";
  
  if ($status) {
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return (defined $sql) ? $self->generic_fetch($sql) : [];
  
}

=head2 fetch_all_by_product_FeatureSet_type

  Arg [1]    : string - product feature set type for this data_set e.g. 'annotated', 'regulatory'
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my $dsets = $dset_adaptor->fetch_all_by_product_FeatureSet_type('regulatory');
  Description: Fetch all DataSets of a given product feature set type
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if no product feaure set type passed
  Caller     : General
  Status     : At Risk - not yet implmented

=cut

sub fetch_all_by_product_FeatureSet_type {
  my ($self, $type, $status) = @_;
  
  throw('Not yet implemented');
  #left join? same for result_set?
  #Or do we do a sneaky workaround and use the FeatureSet adaptor here to get the feature sets we require
  #and then call self?
  #or do we just use feature_sets for reg feats and forget about the possiblity of displaying the suppoorting sets?
  #will we ever want to display supporting sets in contig/neighbourhood view?
  #or will we only use expanded view in a feature context, hence we can just use the Regulatory Attribs?
  #don't really want to expand regulatory feats on contig view, so access FeatureSet directly

  throw("Must provide a supporting_set type argument") if (! defined $type);
  
  my $sql = "ss.type='".$type."'";
  
  if ($status) {
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return (defined $sql) ? $self->generic_fetch($sql) : [];
  
}


=head2 fetch_by_product_FeatureSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : my @dsets = $fs_adaptopr->fetch_by_product_FeatureSet($fset);
  Description: Retrieves DataSet objects from the database based on the FeatureSet.
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : Throws if arg is not a valid FeatureSet
  Caller     : General
  Status     : Deprecated - use fetch_all_by_product_FeatureSet

=cut

sub fetch_all_by_FeatureSet {
  my $self = shift;

  deprecate('Use fetch_by_product_FeatureSet');

  return $self->fetch_by_product_FeatureSet(@_);

}


=head2 fetch_by_product_FeatureSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : my @dsets = $fs_adaptopr->fetch_by_product_FeatureSet($fset);
  Description: Retrieves DataSet objects from the database based on the FeatureSet.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if arg is not a valid FeatureSet
  Caller     : General
  Status     : At Risk

=cut


#This is main FeatureSet, i.e. the result of the analysis of the supporting_sets
#Supporting sets could also be FeatureSets!!!  Confusion!

sub fetch_by_product_FeatureSet {
  my $self = shift;
  my $fset = shift;

    
  if (! ($fset && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet") && $fset->dbID())) {
	throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::FeatureSet object");
  }
	
  my $sql = "ds.feature_set_id = '".$fset->dbID()."'";


  return $self->generic_fetch($sql)->[0];	
}


=head2 fetch_all_by_ResultSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::ResultSet
  Example    : my @dsets = $fs_adaptopr->fetch_all_by_ResultSet($rset);
  Description: Retrieves DataSet objects from the database based on the ResultSet.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if arg is not a valid ResultSet
  Caller     : General
  Status     : At Risk - to be removed

=cut

sub fetch_all_by_ResultSet {
  my $self = shift;
  my $rset = shift;

  deprecate('Use fetch_all_by_supporting_set');

  return $self->fetch_all_by_supporting_set($rset);


  #if(! ($rset && $rset->isa("Bio::EnsEMBL::Funcgen::ResultSet") && $rset->dbID())){
  #  throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::ResultSet object");
  #}
	

  ##self join here to make sure we get all linked result_sets
  #my $sql = 'ds.data_set_id IN (SELECT ds.data_set_id from data_set ds where result_set_id='.$rset->dbID().')';


  #return $self->generic_fetch($sql);	
}



=head2 fetch_all_by_supporting_set

  Arg [1]    : Bio::EnsEMBL::Funcgen::Result|FeatureSet
  Example    : my @dsets = $fs_adaptopr->fetch_all_by_supporting_set($rset);
  Description: Retrieves DataSet objects from the database based on the
               given supporting Result or FeatureSet.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if arg is not a valid Result|FeatureSet
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_supporting_set {
  my $self = shift;
  my $set = shift;

  if (! (ref($set) && 
		 ( $set->isa("Bio::EnsEMBL::Funcgen::ResultSet") || 
		   $set->isa("Bio::EnsEMBL::Funcgen::FeatureSet") ||
		   $set->isa("Bio::EnsEMBL::Funcgen::InputSet")) 
		 && $set->dbID())) {
	throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::ResultSet, FeatureSet or InputSet object");
  }
	
  #self join here to make sure we get all linked result_sets
  my $sql = ' ds.data_set_id IN (SELECT data_set_id from supporting_set where type="'.$set->set_type.'" and supporting_set_id='.$set->dbID().')';
	
  return $self->generic_fetch($sql);	
}





=head2 fetch_all_by_feature_type_class

  Arg [1]    : string - class of associated FeatureSet FeatureType
  Arg [2]    : optional: string - status e.g. 'DISPLAYABLE'
  Example    : my @dsets = @{$ds_adaptopr->fetch_all_by_feature_type_class('HISTONE')};
  Description: Retrieves DataSet objects from the database based on the FeatureSet FeatureType class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if no class arg defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_feature_type_class {
  my ($self, $class, $status) = @_;
  
  throw ('Must provide a FeatureType class to retrieve DataSets') if ! defined $class;
  
  my ($constraint, @dsets);

  if ($status) {
	$constraint = $self->status_to_constraint($status) if $status;
	return [] if ! defined $constraint;
  }


  #This is fetching all feature sets!
  #we need to left join this?
  #we can't do it for class
  #but we can do it for product feature_set type

  foreach my $dset (@{$self->generic_fetch($constraint)}) {
	#uc both here to avoid case sensitivities
	push @dsets, $dset if uc($dset->product_FeatureSet->feature_type->class()) eq uc($class);
  }

  return \@dsets;	
}

=head2 fetch_all_displayable_by_feature_type_class

  Arg [1]    : string - class of associated FeatureSet FeatureType
  Example    : my @dsets = @{$ds_adaptopr->fetch_all_displayable_by_feature_type_class('HISTONE')};
  Description: Wrapper method, retrieves all displayable DataSets with given FeatureSet FeatureType class
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
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
  Returntype : List
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _tables {
  my $self = shift;
  
  return (
		  [ 'data_set',    'ds' ],
		  [ 'supporting_set', 'ss'],
		 );
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  #my $self = shift;

  #will this work? May have multiple record/result_set_id
  
  
  return qw(
			ds.data_set_id	    ds.feature_set_id
			ds.name             ss.type
			ss.supporting_set_id
		   );	
}


=head2 _left_join

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _left_join {
  my $self = shift;
	
  return (['supporting_set', 'ds.data_set_id = ss.data_set_id']);
}


#=head2 _final_clause
#
#  Args       : None
#  Example    : None
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns an ORDER BY clause. Sorting by oligo_feature_id would be
#			   enough to eliminate duplicates, but sorting by location might
#			   make fetching features on a slice faster.
#  Returntype : String
#  Exceptions : None
#  Caller     : generic_fetch
#  Status     : At Risk
#
#=cut


#this should be another left join? on feature_set and a join on feature_type so we can sort lexically on class, name
#should we implement a default sort in the data_set adaptor which could be superceeded by a custom list?

#sub _final_clause {
#	return ' ORDER BY fs.feature_type_id, fs.cell_type_id'; #group on cell_type_id
#}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;
  
  my (@data_sets, @supporting_sets, $data_set, $dbID, $set_id);
  my ($fset_id, $fset, $set, $name, $ss_type, $ss_id);
  
  my %set_adaptors = (
					  feature      => $self->db->get_FeatureSetAdaptor(),
					  result       => $self->db->get_ResultSetAdaptor(),
					  input        => $self->db->get_InputSetAdaptor(),
					 );

  $sth->bind_columns(\$dbID, \$fset_id, \$name, \$ss_type, \$ss_id);
  
  while ( $sth->fetch() ) {
	
    if ((! $data_set) || ($data_set->dbID() != $dbID)) {

	  if ($data_set) {
		$data_set->add_supporting_sets(\@supporting_sets);
		push @data_sets, $data_set;
		#do not set to empty array as this will cause failure of check in DataSet->new
		undef @supporting_sets;
	  }

	  #handle absent sets, dbIDs of 0
      $fset = ($fset_id) ? $set_adaptors{'feature'}->fetch_by_dbID($fset_id) : undef;

      $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(
													  -DBID                => $dbID,
													  -NAME                => $name,
													  -FEATURE_SET         => $fset,
													  -ADAPTOR             => $self,
													 );

	}
	
	if ($ss_id) {	  
	  my $sset = $set_adaptors{$ss_type}->fetch_by_dbID($ss_id);
	  
	  if (! $sset) {
		throw("Could not find supporting $ss_type set with dbID $ss_id");
	  }
	  
	  push @supporting_sets, $sset;
	}
  }

  
  #handle last set
  if ($data_set) {
	#	warn "ssets are @supporting_sets";
	$data_set->add_supporting_sets(\@supporting_sets);
	push @data_sets, $data_set;
  }


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

  Args       : List of Bio::EnsEMBL::Funcgen::DataSet objects
  Example    : $dsa->store(@dsets);
  Description: Stores given DataSet objects in the database. Sets dbID and adaptor 
               on the objects that it stores.
  Returntype : ARRAYREF of stored DataSet objects
  Exceptions : Throws if no DataSet objects passed
               Throws if DataSet object has already been stored
               Throws if any supporting sets have not been stored
  Caller     : General
  Status     : At Risk

=cut

  sub store{
	my ($self, @dsets) = @_;

	throw('Must pass a list of DataSet objects to store') if(! @dsets || $#dsets < 0);


	my $sth = $self->prepare("INSERT INTO data_set (feature_set_id, name) 
                            VALUES (?, ?)");
	my $sth2 = $self->prepare("INSERT INTO supporting_set (data_set_id, supporting_set_id, type) 
                            VALUES (?, ?, ?)");

	my ($fset_id);

	my $db = $self->db();


	foreach my $dset (@dsets) {

	  throw('Must pass a DataSet object to store') if( ! ( ref $dset && 
														   $dset->isa('Bio::EnsEMBL::Funcgen::DataSet')));

	  if ( $dset->is_stored($db) ) {
		throw('DataSet [' . $dset->dbID() . '] is already stored in the database,'.
			  'use store_updated_sets method to add new supporting sets in this DataSet');
	  }
		
	  $fset_id = (defined $dset->product_FeatureSet()) ? $dset->product_FeatureSet->dbID() : 0;
	
	  $sth->bind_param(1, $fset_id,                     SQL_INTEGER);
	  $sth->bind_param(2, $dset->name(),                SQL_VARCHAR);
	  $sth->execute();
	  $dset->dbID( $sth->{'mysql_insertid'} );
	  $dset->adaptor($self);
	

    
	  foreach my $sset (@{$dset->get_supporting_sets()}) {
		
		throw("All supporting Feature and ResultSets must be stored previously.".
			  " Use store_updated_sets method if your DataSet has been stored") if(! $sset->is_stored($db));

		$sth2->bind_param(1, $dset->dbID(),                SQL_INTEGER);
		$sth2->bind_param(2, $sset->dbID(),                SQL_INTEGER);
		$sth2->bind_param(3, $sset->set_type(),            SQL_VARCHAR); #enum feature/result/experimental
		$sth2->execute();
	  }
	}
      
	return \@dsets
  }


=head2 store_updated_sets

  Args       : List of previously stored Bio::EnsEMBL::Funcgen::DataSet objects
  Example    : $dsa->store_updated_sets(@dsets);
  Description: Updates added supporting sets for a given previously stored DataSet
  Returntype : ARRAYREF of updated DataSet objects
  Exceptions : Throws if a list of DataSet objects is not provided
               Throws if DataSet has not been previosuly stored
               Throws if supporting set has not been previously stored
               ? should we throw or warn if a set has been deleted?
  Caller     : General
  Status     : At Risk

=cut

  #This needs to cahnge to an arrayref of dset and an overwrite flag


  sub store_updated_sets{
	my ($self, $dsets, $overwrite) = @_;

	throw('Must pass a list of DataSet objects to store') if(! @$dsets || $#{$dsets} < 0);
	my ($sql);
	my $db = $self->db();
	my $sth = $self->prepare("INSERT INTO supporting_set (data_set_id, supporting_set_id, type) 
                            VALUES (?, ?, ?)");

	foreach my $dset (@$dsets) {
	  throw('Must pass a DataSet object to update') if( ! ( ref $dset && 
															$dset->isa('Bio::EnsEMBL::Funcgen::DataSet')));
	
	  throw('DataSet [' . $dset->dbID() . '] must be previous stored in the database') if (! $dset->is_stored($db) );
	  my $stored_dset = $self->fetch_by_name($dset->name);


	  #Update product FeatureSet
	  #We need to do this first so we cacn check wether we're updated supporting_sets
	  #for a data set which has already got a product FeatureSet...not wise

	  my $fset = $dset->product_FeatureSet;
	  my $stored_fset = $stored_dset->product_FeatureSet;
	  #This fset check is slight overkill, as you have to severly mangle a dataset to fail this validation

	  if (defined $stored_fset) {

		if (! defined $fset) {
		  #How will this have ever happened?
		  warn("Populating absent product FeatureSet from DB for DataSet:\t".$dset->name);
		} else {
		  #validate sets
		  if ($fset->dbID != $stored_fset->dbID) {
			my $msg = 'Found product FeatureSet mismatch whilst updating DataSet('.$dset->name.
			  "):\tStored:".$stored_fset->name."\tUpdate:".$fset->name;
			throw($msg) if ! $overwrite;
			warn $msg;
		  }
		}
	  } else {
		#update data_set table
		$sql = 'update data_set set feature_set_id='.$fset->dbID.' where data_set_id='.$dset->dbID;
		$self->dbc->do($sql);
	  }

	  my @sorted_ssets = sort {$a->dbID <=> $b->dbID} @{$dset->get_supporting_sets};
	  my @stored_ssets = sort {$a->dbID <=> $b->dbID} @{$stored_dset->get_supporting_sets};
	  my $mismatch = 0;
	
	  $mismatch = 1 if(scalar(@sorted_ssets) != scalar(@stored_ssets));
	
	  if (! $mismatch) {
	  
		for my $i (0..$#stored_ssets) {
		
		  if ($stored_ssets[$i]->dbID != $sorted_ssets[$i]->dbID) {
			$mismatch=1;
			last;
		  }
		}
	  }

	  #Delete old supporting_sets
	  #We really only want to do this if there is a difference
	  #batched jobs cause race condition here
	  #unless updated once before submission	
	  if ($mismatch &&
		  $overwrite) {
		$sql = 'DELETE from supporting_set where data_set_id='.$dset->dbID;
		eval { $self->db->dbc->do($sql) };
		throw("Couldn\'t delete supporting_sets for DataSet:\t".$stored_dset->name."\n$@") if $@;
		@stored_ssets = ();
	  }


	  #Update supporting sets
	  my %stored_dbids;
	  map {$stored_dbids{$_->dbID} = undef} @stored_ssets if @stored_ssets;

	  foreach my $sset (@sorted_ssets) {	
		my $dbid = $sset->dbID;

		if (! grep(/^${dbid}$/, keys %stored_dbids)) {		
		  throw("All supporting sets must be stored previously.") if(! $sset->is_stored($db));

		  #This causes problems when we want to re-run by slice
		  #Currently need
		  #warn "temporarily suspended throw for existing feature_set";
		  throw('You are trying to update supporting sets for a data set which already has a product FeatureSet('.$stored_fset->name.').  Either rollback the FeatureSet before adding more supporting sets or specify the overwrite flag.') if defined $stored_fset && ! $overwrite;

		  $sth->bind_param(1, $dset->dbID,            SQL_INTEGER);
		  $sth->bind_param(2, $sset->dbID,            SQL_INTEGER);
		  $sth->bind_param(3, $sset->set_type,        SQL_VARCHAR);
		  #How is this failing?
		  #Is the index not being updated after the delete
		  $sth->execute();
		}
	  }
	}
      
	return $dsets
  }


  1;

