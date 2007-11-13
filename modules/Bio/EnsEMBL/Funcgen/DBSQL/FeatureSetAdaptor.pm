#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor - A database adaptor for fetching and
storing Funcgen feature sets.

=head1 SYNOPSIS

my $fs_adaptor = $db->get_FeatureSetAdaptor();

my @fsets = $fs_adaptor->fetch_all_by_Experiment($exp);
my @displayable_fsets = @{$fs_adaptor->fetch_all_displayable()};

=head1 DESCRIPTION

The FeatureSetAdaptor is a database adaptor for storing and retrieving
Funcgen feature set.  

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);



=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_FeatureType($type);
  Description: Retrieves FeatureSet objects from the database based on feature_type id.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if arg is not a valid FeatureType
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType {
    my $self = shift;
    my $ftype = shift;
    my $status = shift;
    
    if(! ($ftype && $ftype->isa("Bio::EnsEMBL::Funcgen::FeatureType"))){
      throw("Must provide a valid Bio::EnsEMBL::Funcgen::FeatureType object");
    }
	
    my $sql = "fs.feature_type_id = '".$ftype->dbID()."'";

    if($status){
      my $constraint = $self->status_to_constraint($status);
      $sql = (defined $constraint) ? $sql." ".$constraint : undef;
    }

    return $self->generic_fetch($sql);	
}


=head2 fetch_all_by_type

  Arg [1]    : String - Type of feature set i.e. 'annotated', 'regulatory' or 'supporting'
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_type('annotated');
  Description: Retrieves FeatureSet objects from the database based on feature_set type.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if type not defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_type {
    my $self = shift;
    my $type = shift;
    my $status = shift;
    
    throw('Must provide a feature_set type') if(! defined $type);
	
    my $sql = "fs.type = '".$type."'";

    if($status){
      my $constraint = $self->status_to_constraint($status);
      $sql = (defined $constraint) ? $sql." AND ".$constraint : undef;
    }

    return $self->generic_fetch($sql);	
}

=head2 fetch_all_displayable_by_type

  Arg [1]    : String - Type of feature set i.e. 'annotated', 'regulatory' or 'supporting'
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_type('annotated');
  Description: Wrapper method for fetch_all_by_type
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_displayable_by_type {
    my $self = shift;
    my $type = shift;
  
	$self->fetch_all_by_type($type, 'DISPLAYABLE');
}


=head2 fetch_all_by_CellType

  Arg [1]    : Bio::EnsEMBL::Funcgen::CellType
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_CellType($ctype);
  Description: Retrieves FeatureSet objects from the database based on the CellType.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if arg is not a valid CellType
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_CellType {
    my $self = shift;
    my $ctype = shift;
    my $status = shift;
    
    if(! ($ctype && $ctype->isa("Bio::EnsEMBL::Funcgen::CellType") && $ctype->dbID())){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::CellType object");
    }
	
    my $sql = "fs.cell_type_id = '".$ctype->dbID()."'";

    if($status){
      my $constraint = $self->status_to_constraint($status) if $status;
      $sql = (defined $constraint) ? $sql." ".$constraint : undef;
    }

    return $self->generic_fetch($sql);	
}



=head2 fetch_all_by_FeatureType_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : Bio::EnsEMBL::Analysis
  Arg [3]    : (optional) Bio::EnsEMBL::Funcgen::CellType
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_FeatureType_Analysis($ftype, $anal, $ctype);
  Description: Retrieves FeatureSet objects from the database based on FeatureType, Analysis and 
               CellType if defined.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if args 1 and 2 are not valid or stored
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType_Analysis {
    my ($self, $ftype, $anal, $ctype) = @_;
    

	my $sql = '';

    if(! ($ftype && $ftype->isa("Bio::EnsEMBL::Funcgen::FeatureType") && $ftype->dbID())){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::FeatureType object");
    }
	
	if(! ($anal && $anal->isa("Bio::EnsEMBL::Analysis") && $anal->dbID())){
      throw("Must provide a valid stored Bio::EnsEMBL::Analysis object");
    }

	if($ctype){

	  if(! ($ctype->isa("Bio::EnsEMBL::Funcgen::CellType") && $ctype->dbID())){
		throw("Argument must be a valid stored Bio::EnsEMBL::Funcgen::CellType object");
	  }
	  
	  $sql = ' AND fs.cell_type_id='.$ctype->dbID();
	}


    $sql = 'fs.feature_type_id ='.$ftype->dbID().' AND fs.analysis_id='.$anal->dbID().$sql;

  
    return $self->generic_fetch($sql);	
}

=head2 fetch_by_name

  Arg [1]    : string - name of FeatureSet
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = @{$fset_adaptor->fetch_by_name('feature_set-1')};
  Description: Fetch all FeatureSets wit a given name
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if no name passed 
  Caller     : General
  Status     : At Risk - change to fetch_by_name when name is made unique key

=cut

sub fetch_by_name {
  my ($self, $name, $status) = @_;
  
  throw("Must provide a name argument") if (! defined $name);
  
  my $sql = "fs.name='".$name."'";
  
  if($status){
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return $self->generic_fetch($sql)->[0];
  
}

#=head2 fetch_by_external_db

#  Arg [1]    : string - name of external_db
#  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
#  Example    : my @fsets = @{$fset_adaptor->fetch_by_external_db('miRanda')};
#  Description: Fetch all FeatureSets which are linked to an external_db of a given name
#  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet objects
#  Exceptions : Throws if no external_db name passed 
#  Caller     : General
#  Status     : At Risk 

#=cut

#sub fetch_by_external_db {
#  my ($self, $name, $status) = @_;
  
#  throw("Must provide a name argument") if (! defined $name);
  
#  my $sql = "ed.name='".$name."'";
  
#  if($status){
#    my $constraint = $self->status_to_constraint($status) if $status;
#    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
#  }

#  return $self->generic_fetch($sql)->[0];
  
#}


=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoArray - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty OligoArray objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::OligoArray getters
  Status     : Medium Risk

=cut

sub fetch_attributes {
    my $self = shift;
    my $array = shift;

    my $tmp_array = $self->fetch_by_dbID( $array->dbID() );
    %$array = %$tmp_array;
}

=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _tables {
	my $self = shift;
	
	return (['feature_set',     'fs'],
			#['feature_set_db', 'fsd'],
			#['external_db',     'ed'],
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
	
	return qw( fs.feature_set_id fs.feature_type_id fs.analysis_id fs.cell_type_id fs.name fs.type);# ed.db_name);
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

#sub _left_join {
#  my $self = shift;
	
#  return (['feature_set_db', 'fs.feature_set_id = fsd.feature_set_id']);
#}

=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : String
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

#sub _default_where_clause {
#  my $self = shift;
  #will this return if there are no entrie in data_set_member?
  #do we have to implement a join here?

	
  #return 'fsd.external_db_id = ed.external_db_id';
#}



=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@fsets, $fset, $analysis, %analysis_hash, $feature_type, $cell_type, $name, $type);#, $dbname);
	my ($feature_set_id, $ftype_id, $analysis_id, $ctype_id, %ftype_hash, %ctype_hash);
	
	my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
	my $anal_adaptor = $self->db->get_AnalysisAdaptor();
	my $ct_adaptor = $self->db->get_CellTypeAdaptor();
	$ctype_hash{'NULL'} = undef;

	$sth->bind_columns(\$feature_set_id, \$ftype_id, \$analysis_id, \$ctype_id, \$name, \$type);# \$dbname);
	
	while ( $sth->fetch()) {

		$ctype_id ||= 'NULL';

		# Get the analysis object
		$analysis_hash{$analysis_id} = $anal_adaptor->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});

		# Get the feature type object
		$ftype_hash{$ftype_id} = $ft_adaptor->fetch_by_dbID($ftype_id) if(! exists $ftype_hash{$ftype_id});
		
		# Get the cell_type object
		$ctype_hash{$ctype_id} = $ct_adaptor->fetch_by_dbID($ctype_id) if(! exists $ctype_hash{$ctype_id});


		$fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
													   -dbID         => $feature_set_id,
													   -adaptor      => $self,
													   -feature_type => $ftype_hash{$ftype_id},
													   -analysis     => $analysis_hash{$analysis_id},
													   -cell_type    => $ctype_hash{$ctype_id},
													   -name         => $name,
													   -type         => $type,
													   #-external_db_name  => $dbname,
							      );

		push @fsets, $fset;

	}

	return \@fsets;
}




=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Example    : $oaa->store($fset1, $fset2, $fset3);
  Description: Stores FeatureSet objects in the database.
  Returntype : Listref of stored FeatureSet objects
  Exceptions : Throws if FeatureSet does not have a stored FeatureType
               Throws if invalid FeatureSet passed
               Warns if external_db_name not defined is type is external
               Throws if external_db is not present in the db
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @fsets = @_;

	#if(scalar(@fsets) == 0);

	my $sth = $self->prepare("INSERT INTO feature_set
                                 (feature_type_id, analysis_id, cell_type_id, name, type)
                                 VALUES (?, ?, ?, ?, ?)");


	my $esd_sth = $self->prepare("INSERT INTO feature_set_db
                                 (feature_set_id, external_db_id)
                                 VALUES (?, ?)");

	my ($sql, $edb_id, %edb_hash);

    foreach my $fset (@fsets) {
		throw('Can only store FeatureSet objects, skipping $fset')	if ( ! $fset->isa('Bio::EnsEMBL::Funcgen::FeatureSet'));
		
		if (!( $fset->dbID() && $fset->adaptor() == $self )){#use is_stored?

		  #this should use the xref API directly
		  
		  #if($fset->type() eq 'external' && ! defined $fset->external_db_name()){
		  #	warn('You are loading an ExternalFeature FeatureSet with no associated external_db name');
		  #  }


		  #Need to check external_db is present.
		#  if(defined $fset->external_db_name() && ! exists $edb_hash{$fset->external_db_name()}){
		#	$sql = 'SELECT external_db_id from external_db where db_name="'.$fset->external_db_name().'"';
		#	($edb_id) = $self->db->dbc->db_handle->selectrow_array($sql);
		##	
		#	throw ('You must specifcy a previously stored external_db name') if(! $edb_id);
		#	$edb_hash{$fset->external_db_name()} = $edb_id;
		  #}
					
			#my $s_fset = $self->fetch_by_unique_and_experiment_id($ec->unique_id(), $ec->experiment_id());
			#throw("ExperimentalChip already exists in the database with dbID:".$s_ec->dbID().
			#	  "\nTo reuse/update this ExperimentalChip you must retrieve it using the ExperimentalChipAdaptor".
			#	  "\nMaybe you want to use the -recover option?") if $s_ec;
					 
			throw("FeatureSet must have a stored FeatureType") if (! $fset->feature_type->is_stored($self->db()));
			 
			my $ctype_id = (defined $fset->cell_type()) ? $fset->cell_type->dbID() : undef;

			$sth->bind_param(1, $fset->feature_type->dbID(), SQL_INTEGER);
			$sth->bind_param(2, $fset->analysis->dbID(),     SQL_INTEGER);
			$sth->bind_param(3, $ctype_id,                   SQL_INTEGER);
			$sth->bind_param(4, $fset->name(),               SQL_VARCHAR);
			$sth->bind_param(5, $fset->type(),               SQL_VARCHAR);
		
			$sth->execute();
			$fset->dbID($sth->{'mysql_insertid'});
			$fset->adaptor($self);

			#if(defined $edb_id){
			#  $esd_sth->bind_param(1, $fset->dbID(), SQL_INTEGER);
			#  $esd_sth->bind_param(2, $edb_id,       SQL_INTEGER);
			#  $esd_sth->execute();
			#}

		}else{
			#assume we want to update the states
			warn('You may want to use $fset->adaptor->store_states($fset)');
			$self->store_states($fset);
		}
	}
	return \@fsets;
}

=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$oaa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoArray objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('feature_set');
}



1;

