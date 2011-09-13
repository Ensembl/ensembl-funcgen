#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor
#

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

Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor - A database adaptor for fetching and
storing Funcgen feature sets.

=head1 SYNOPSIS

my $fs_adaptor = $db->get_FeatureSetAdaptor();

my @fsets = $fs_adaptor->fetch_all_by_Experiment($exp);
my @displayable_fsets = @{$fs_adaptor->fetch_all_displayable()};

=head1 DESCRIPTION

The FeatureSetAdaptor is a database adaptor for storing and retrieving
Funcgen feature set.  

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
      $sql = (defined $constraint) ? $sql." AND ".$constraint : undef;
    }

    return (defined $sql) ? $self->generic_fetch($sql) : [];	
}


=head2 fetch_all_by_type

  Arg [1]    : String - Type of feature set i.e. 'annotated', 'regulatory' or 'external'
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
    
	
	
    return $self->fetch_all_by_feature_class($type, $status);	
}

=head2 fetch_all_by_feature_class

  Arg [1]    : String - feature class i.e. 'annotated', 'regulatory' or 'external'
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_feature_class('annotated');
  Description: Retrieves FeatureSet objects from the database based on feature_set type.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if type not defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_feature_class {
    my $self = shift;
    my $type = shift;
    my $status = shift;
    
    throw('Must provide a feature_set type') if(! defined $type);
	
    my $sql = "fs.type = '".$type."'";

    if($status){
      my $constraint = $self->status_to_constraint($status);
      $sql = (defined $constraint) ? $sql." AND ".$constraint : undef;
    }

    return (defined $sql) ? $self->generic_fetch($sql) : [];	
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

    return (defined $sql) ? $self->generic_fetch($sql) : [];
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

  return (defined $sql) ? $self->generic_fetch($sql)->[0] : [];
  
}


=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureSet - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty FeatureSet objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::FeatureSet getters
  Status     : At Risk - Not implemented for FeatureSets, remove?

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
	
	return (['feature_set',     'fs']);
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
	
	return qw( fs.feature_set_id fs.feature_type_id fs.analysis_id fs.cell_type_id fs.name fs.type fs.description fs.display_label fs.experiment_id);
}




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
	
	my (@fsets, $fset, $analysis, %analysis_hash, $feature_type, $cell_type, $name, $type, $display_label, $desc);
	my ($feature_set_id, $ftype_id, $analysis_id, $ctype_id, $exp_id, %ftype_hash, %ctype_hash);
	
	my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
	my $anal_adaptor = $self->db->get_AnalysisAdaptor();
	my $ct_adaptor = $self->db->get_CellTypeAdaptor();
	$ctype_hash{'NULL'} = undef;

	$sth->bind_columns(\$feature_set_id, \$ftype_id, \$analysis_id, \$ctype_id, \$name, \$type, \$desc, \$display_label, \$exp_id);
	
	while ( $sth->fetch()) {

		$ctype_id ||= 'NULL';

		# Get the analysis object
		$analysis_hash{$analysis_id} = $anal_adaptor->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});

		# Get the feature type object
		$ftype_hash{$ftype_id} = $ft_adaptor->fetch_by_dbID($ftype_id) if(! exists $ftype_hash{$ftype_id});
		
		# Get the cell_type object
		$ctype_hash{$ctype_id} = $ct_adaptor->fetch_by_dbID($ctype_id) if(! exists $ctype_hash{$ctype_id});
		
		#Use new_fast here and strip the prefixed -'s
		$fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
		  (
		   -dbID          => $feature_set_id,
		   -adaptor       => $self,
		   -feature_type  => $ftype_hash{$ftype_id},
		   -analysis      => $analysis_hash{$analysis_id},
		   -cell_type     => $ctype_hash{$ctype_id},
		   -name          => $name,
		   -feature_class => $type,
		   -display_label => $display_label,
		   -description   => $desc,
		   -experiment_id => $exp_id,
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
               Throws if not FeatureSets passed
               Warns if external_db_name not defined is type is external
               Throws if external_db is not present in the db
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @fsets = @_;

	throw('Must supply a list of FeatureSets to store') if(scalar(@fsets) == 0);

	my $sth = $self->prepare("INSERT INTO feature_set
                                 (feature_type_id, analysis_id, cell_type_id, name, type, description, display_label, experiment_id)
                                 VALUES (?, ?, ?, ?, ?, ?, ?, ?)");

	my $db = $self->db;
	my ($sql, $edb_id, %edb_hash);
	
    foreach my $fset (@fsets) {
		throw('Can only store FeatureSet objects, skipping $fset')	if ( ! $fset->isa('Bio::EnsEMBL::Funcgen::FeatureSet'));
		
		if (! $fset->is_stored($db) ) {

		  # Check FeatureType and Analysis
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fset->feature_type);
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $fset->analysis);
			 

		  # Check optional Experiment and CellType
		  my $ctype_id;
		  my $ctype = $fset->cell_type;

		  if($ctype){
			$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $ctype);
			$ctype_id = $ctype->dbID;
		  }

		  my $exp_id;
		  my $exp =  $fset->get_Experiment; 

		  if($exp){
			$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exp);
			$exp_id = $exp->dbID;
		  }
		  
		  $sth->bind_param(1, $fset->feature_type->dbID,     SQL_INTEGER);
		  $sth->bind_param(2, $fset->analysis->dbID,         SQL_INTEGER);
		  $sth->bind_param(3, $ctype_id,                     SQL_INTEGER);
		  $sth->bind_param(4, $fset->name,                   SQL_VARCHAR);
		  $sth->bind_param(5, $fset->feature_class,          SQL_VARCHAR);
		  $sth->bind_param(6, $fset->description,            SQL_VARCHAR);
		  $sth->bind_param(7, $fset->display_label,          SQL_VARCHAR);
		  $sth->bind_param(8, $exp_id,                       SQL_INTEGER);
		  
		  		  
		  $sth->execute();
		  $fset->dbID($sth->{'mysql_insertid'});
		  $fset->adaptor($self);
		}
		else{
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


=head2 fetch_focus_set_config_by_FeatureSet

  Args       : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : $self->{'focus_set'} = $self->adaptor->fetch_focus_set_config_by_FeatureSet($self);
  Description: Caches and returns focus set config for a given FeatureSet
  Returntype : Boolean
  Exceptions : Warns if meta entry not present
  Caller     : Bio::EnsEMBL::Funcgen::FeatureSet::is_focus_set
  Status     : At Risk

=cut

sub fetch_focus_set_config_by_FeatureSet{
  my ($self, $fset) = @_;
  
  $self->{focus_set_config} ||= {};
  
  if(! defined $self->{focus_set_config}->{$fset->dbID}){
	
	#Is is an attribute set?
	if($self->fetch_attribute_set_config_by_FeatureSet($fset)){
	  
	  #Need to define these as RegBuild config
	  if( ($fset->feature_type->class eq 'Transcription Factor') ||
		  ($fset->feature_type->class eq 'Open Chromatin') ){
		$self->{focus_set_config}->{$fset->dbID} = 1;
	  }
	}
  }
	
  return $self->{focus_set_config}->{$fset->dbID};
}


=head2 fetch_attribute_set_config_by_FeatureSet

  Args       : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : $self->{'attribute_set'} = $self->adaptor->fetch_attribute_set_config_by_FeatureSet($self);
  Description: Caches and returns attribute set config for a given FeatureSet
  Returntype : Boolean
  Exceptions : Warns if meta entry not present
  Caller     : Bio::EnsEMBL::Funcgen::FeatureSet::is_attribute_set
  Status     : At Risk

=cut

sub fetch_attribute_set_config_by_FeatureSet{
    my ($self, $fset) = @_;

	$self->{attribute_set_config} ||= {};

	if(! defined $self->{attribute_set_config}->{$fset->dbID}){
	  $self->{attribute_set_config}->{$fset->dbID} = 0;  #set cache default
	  my $meta_key =  'regbuild.'.$fset->cell_type->name.'.feature_set_ids';

	  #list_value_by_key caches, so we don't need to implement this in the adaptor
	  my ($attr_ids) = @{$self->db->get_MetaContainer->list_value_by_key($meta_key)};

	  if(! defined $attr_ids){
		warn("Cannot detect attribute set as meta table does not contain $meta_key");
	  }
	  else{

		foreach my $aid(split/,\s*/, $attr_ids){
		  $self->{attribute_set_config}->{$aid} = 1;
		}
	  }
	}

    return $self->{attribute_set_config}->{$fset->dbID};
  }




1;

