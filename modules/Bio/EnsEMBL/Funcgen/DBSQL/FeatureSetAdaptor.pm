#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor - A database adaptor for fetching and
storing Funcgen feature sets.

=head1 SYNOPSIS

my $fs_adaptor = $db->get_FeatureSetAdaptor();


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
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);
#use base does not import %true_tables or %tables or in fact %sql_types
#so cannot use base for any of the adaptors

#Exported from BaseAdaptor
#$true_tables{feature_set} = [  [ 'feature_set', 'fs' ] ];
#This derefs as we push onto tables
#and we don't want to alter true table
#@{$tables{feature_set}} = @{$true_tables{feature_set}};


#Look into using ReadOnly::array?
#This will provide true read only, but will need to be our'd and export/imported
#use in conjuction?
#use Readonly;
#Readonly::Array my @true_tables => ([ 'feature_set', 'fs' ]);
#use constant TRUE_TABLES => [ @true_tables ];
#does not work

use constant TRUE_TABLES => [[ 'feature_set', 'fs' ]];
use constant TABLES      => [[ 'feature_set', 'fs' ]];

#Currently these need to be listrefs [[], ...] and push directly onto TABLE rather than _tables




#@{$tables{feature_set}} = (TRUE_TABLES);
#use constant TABLES => (TRUE_TABLES); #this does not deref, hence pushes affect TRUE_TABLES!


#use constant here still allows the contents of the ref to be modified
#Simply prevents need for import/export

#Now change %tables to an attribute!
#Can't set as an attribute here, would have to be done in new or _tables


#Had to use hashes to prevent different adaptors resetting package level global vars
#TRUE_TABLES does not require this.


#No need for true_final_clause
	

=head2 fetch_all

  Arg [1]    : optional HASHREF - Parameter hash containing contraints config lists e.g.
                  {'constraints' => 
                    {
                     cell_types     => [$cell_type, ...],     # Bio::EnsEMBL::Funcgen::CellType
                     feature_types  => [$ftype, ...],         # Bio::EnsEMBL::Funcgen::FeatureType
                     evidence_types => [$evidence_type, ...], # String e.g. 'DNase1 & TFBS' or 'Hists & Pols'
                     status         => $status_name,          # String e.g. IMPORTED
                     analyses       => [$analysis, ...]       # Bio::EnsEMBL::Analysis
                     projects       => [$proj_name, ...]      # String (experimental_group.name is_project=1) e.g. ENCODE
                    }
                  } 
  Example    : 
  Description: Retrieves a list of FeatureSets. Optional parameters hash allows for flexible query terms.
  Returntype : ARRAYREF of Bio::EnsEMBL::FeatureSet objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all{
  my ($self, $params_hash) = @_;

  my $results = $self->generic_fetch($self->compose_constraint_query($params_hash));
  #@{$tables{feature_set}} = @{$true_tables{feature_set}}; #in case we have added tables e.g. status
  $self->reset_true_tables;

  return $results;
}



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
    my ($self, $ftype, $status) = @_;
    
	my $params = {constraints => {feature_types => [$ftype]}};
	$params->{constraints}{status} = $status if $status;
	#No need to reset tables for these
	return $self->generic_fetch($self->compose_constraint_query($params));	
}


=head2 fetch_all_by_type

  Arg [1]    : String - Type of feature set i.e. 'annotated', 'regulatory', 'segmentation' or 'external'
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
    
	#deprecate this?
	
    return $self->fetch_all_by_feature_class($type, $status);	
}


=head2 fetch_all_by_feature_class

  Arg [1]    : String - feature class i.e. 'annotated', 'regulatory', 'segmentation' or 'external'
  Arg [2]    : String (optional) - status e.g. 'DISPLAYABLE'
  Arg [2]    : Bio::EnsEMBL::Funcgen::CellType (optional) or a HASH parameters 
               containing contraint config e.g.

                   $feature_set_adaptor->fetch_all_by_feature_class
                                           ('annotated', 
                                             {'constraints' => 
                                               {
                                               cell_types     => [$cell_type], #Bio::EnsEMBL::Funcgen::CellType
                                               projects       => ['ENCODE'],
                                               evidence_types => ['Hists & Pols'],
                                               feature_types  => [$ftype], #Bio::EnsEMBL::Funcgen::FeatureType
                                               } 
                                             });

  Example    : my @fsets = @{$fs_adaptopr->fetch_all_by_feature_class('annotated')};
  Description: Retrieves FeatureSet objects from the database based on feature_set type.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if type not defined
  Caller     : General
  Status     : At Risk - Move status to params hash

=cut

sub fetch_all_by_feature_class {
  my ($self, $type, $status, $params) = @_;
  
  throw('Must provide a feature_set type') if(! defined $type);	
  my $sql = "fs.type = '".$type."'";
  
  if (defined $params) {        #Some redundancy over $ctype arg and $params cell_type

    if ( ref($params) eq 'Bio::EnsEMBL::Funcgen::CellType') {
      $params = {constraints => {cell_types => [$params]}};
    } elsif (ref($params) ne 'HASH') {
      throw('Argument must be a Bio::EnsEMBL::Funcgen::CellType or a params HASH');
    }
  }


  if ($status) {
    $params->{constraints}{status} = $status;
  }
  
  #Deal with params constraints
  my $constraint = $self->compose_constraint_query($params);
  $sql .=  " AND $constraint " if $constraint;


  #Get result and reset true tables
  my $result = (defined $sql) ? $self->generic_fetch($sql) : [];
  $self->reset_true_tables;

  return $result;
}





=head2 fetch_all_displayable_by_type

  Arg [1]    : String - Type of feature set i.e. 'annotated', 'regulatory' or 'supporting'
  Arg [2]    : Bio::EnsEMBL::Funcgen::CellType (optional) or parameters HASH
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_type('annotated');
  Description: Wrapper method for fetch_all_by_type
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_displayable_by_type {
    my ($self, $type, $ctype_or_params) = @_;
  
	#Move status to config hash
	$self->fetch_all_by_feature_class($type, 'DISPLAYABLE', $ctype_or_params);
	
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
  my ($self, $ctype, $status) = @_;
  
  my $params = {constraints => {cell_types => [$ctype]}};
  $params->{constraints}{status} = $status if ($status);
  my $results = $self->generic_fetch($self->compose_constraint_query($params));
  #@{$tables{feature_set}} = @{$true_tables{feature_set}}; #in case we have added status
  $self->reset_true_tables;
	return $results;	
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
  
  my $params = {constraints => 
				{
				 feature_types => [$ftype],
				 analyses     => [$anal],
				}
			   };

  $params->{constraints}{cell_types} = [$ctype] if $ctype;
  return $self->generic_fetch($self->compose_constraint_query($params));	

}

=head2 fetch_by_name

  Arg [1]    : string - name of FeatureSet
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = @{$fset_adaptor->fetch_by_name('feature_set-1')};
  Description: Fetch all FeatureSets wit a given name
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if no name passed 
  Caller     : General
  Status     : At Risk

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

	#return @{$tables{feature_set}};
  return ( @{$self->TABLES} );
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
	
	return qw( fs.feature_set_id fs.feature_type_id 
			   fs.analysis_id fs.cell_type_id 
			   fs.name fs.type 
			   fs.description fs.display_label 
			   fs.input_set_id);
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
	my ($feature_set_id, $ftype_id, $analysis_id, $ctype_id, $input_set_id, %ftype_hash, %ctype_hash);
	
	my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
	my $anal_adaptor = $self->db->get_AnalysisAdaptor();
	my $ct_adaptor = $self->db->get_CellTypeAdaptor();
	$ctype_hash{'NULL'} = undef;

	$sth->bind_columns(\$feature_set_id, \$ftype_id, \$analysis_id, \$ctype_id, 
                       \$name, \$type, \$desc, \$display_label, \$input_set_id);
	
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
		   -input_set_id  => $input_set_id,
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

	my $sth = $self->prepare
	  (
	   "INSERT INTO feature_set
        (feature_type_id, analysis_id, cell_type_id, name, type, description, display_label, input_set_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
	  );

	my $db = $self->db;
	my ($sql, $edb_id, %edb_hash);
	
    foreach my $fset (@fsets) {
		throw('Can only store FeatureSet objects, skipping $fset')	if ( ! $fset->isa('Bio::EnsEMBL::Funcgen::FeatureSet'));
		
	
		if (! $fset->is_stored($db) ) {

		  # Check FeatureType and Analysis
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fset->feature_type);
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $fset->analysis);
			 

		  # Check optional InputSet and CellType
		  my $ctype_id;
		  my $ctype = $fset->cell_type;

		  if($ctype){
			$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $ctype);
			$ctype_id = $ctype->dbID;
		  }

		  my $input_set_id;
		  my $input_set =  $fset->get_InputSet; 

		  if($input_set){
			$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $input_set);
			$input_set_id = $input_set->dbID;
		  }

		  
		  $sth->bind_param(1, $fset->feature_type->dbID,     SQL_INTEGER);
		  $sth->bind_param(2, $fset->analysis->dbID,         SQL_INTEGER);
		  $sth->bind_param(3, $ctype_id,                     SQL_INTEGER);
		  $sth->bind_param(4, $fset->name,                   SQL_VARCHAR);
		  $sth->bind_param(5, $fset->feature_class,          SQL_VARCHAR);
		  $sth->bind_param(6, $fset->description,            SQL_VARCHAR);
		  $sth->bind_param(7, $fset->display_label,          SQL_VARCHAR);
		  $sth->bind_param(8, $input_set_id,                 SQL_INTEGER);
		  
		  $sth->execute();
		  $fset->dbID($sth->{'mysql_insertid'});
		  $fset->adaptor($self);
		}
		else{
			warn('FeatureSet '.$fset->name.'is already stored, updating status entries');
			$self->store_states($fset);
		}
	}
	return \@fsets;
}

#focus is old terminology
#tend to use 'core' now
#although here we use focus to determine a set which is actually used in the build
#where as core is at the feature type level
#i.e. you can have a feature set with a core feature type 
#which is not in the build and is therefore not a focus set.

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
  
  if (! defined $self->{focus_set_config}->{$fset->dbID}) {
	
    #Is is an attribute set?
    if ($self->fetch_attribute_set_config_by_FeatureSet($fset)) {
	  
      #Need to define these as RegBuild config
      if ($fset->feature_type->is_core_evidence) {
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

	if (! defined $self->{attribute_set_config}->{$fset->dbID}) {
	  $self->{attribute_set_config}->{$fset->dbID} = 0; #set cache default
	  my $string_key =  'regbuild.'.$fset->cell_type->name.'.feature_set_ids';

	  #list_value_by_key caches, so we don't need to implement this in the adaptor
	  #my ($attr_ids) = @{$self->db->get_MetaContainer->list_value_by_key($meta_key)};

	  my $species_id = $self->db->species_id;
	  my ($attr_ids) = $self->db->dbc->db_handle->selectrow_array("SELECT string from regbuild_string where name='${string_key}' and species_id=$species_id");

	  if (! defined $attr_ids) {
      warn("Cannot detect attribute set as regbuild_string table does not contain $string_key");
	  } 
    else {

      foreach my $aid (split/,\s*/, $attr_ids) {
        $self->{attribute_set_config}->{$aid} = 1;
      }
	  }
	}

  return $self->{attribute_set_config}->{$fset->dbID};
}




sub fetch_feature_set_filter_counts{
  my $self = shift;

   my $sql = 'SELECT count(*), eg.name, eg.description, eg.is_project, ft.class, ct.name, ct.description '.
    'FROM experimental_group eg, experiment e, feature_set fs, feature_type ft, cell_type ct, '.
      'status s, status_name sn, input_set inp '.
        'WHERE fs.input_set_id=inp.input_set_id and inp.experiment_id=e.experiment_id '.
          'AND e.experimental_group_id=eg.experimental_group_id '.
          'AND fs.feature_type_id=ft.feature_type_id AND fs.cell_type_id=ct.cell_type_id '.
            'AND fs.feature_set_id=s.table_id AND fs.type="annotated" AND s.table_name="feature_set" '.
              'AND s.status_name_id=sn.status_name_id and sn.name="DISPLAYABLE" '.
                'GROUP BY eg.name, eg.is_project, ft.class, ct.name';

  #warn $sql;
  #Need to write HC around this as we sometimes get less than expect. 
  

  my @rows = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
  my $ft_a = $self->db->get_FeatureTypeAdaptor;
  my $ftype_info = $ft_a->get_regulatory_evidence_info;

  my %filter_info = 
    ( 
     #Project=> {},
     #'Cell/Tissue' => {},
     All =>
     { 
      All =>{ count       => 0,
              description => 'All experiments',
            }
     }
     
    );
  
  foreach my $row(@rows){

    my ($count, $project, $proj_desc, $is_proj, $ft_class, $ct_name, $ct_desc) = @$row;
    
    #All counts
    $filter_info{All}{All}{count} += $count;
  
    #Project counts
    if($is_proj){
      
      if(! exists $filter_info{Project}{$project}){
        $filter_info{Project}{$project} = 
          { count       => 0,
            description => $proj_desc,
          };
      }

      $filter_info{Project}{$project}{count} += $count;
    }

    #Cell/Tissue counts
    if(! exists $filter_info{'Cell/Tissue'}{$ct_name}){
      $filter_info{'Cell/Tissue'}{$ct_name} = 
        { count       => 0,
          description => $ct_desc,
        };
    }   
    $filter_info{'Cell/Tissue'}{$ct_name}{count} += $count;
    
    #Evidence class counts
    #Do we want to split this into ft.class
    #i.e. split 'DNase1 & TFBS'
    my $evidence_type  = $ft_a->get_regulatory_evidence_type($ft_class);
    my $ft_class_label =  $ftype_info->{$evidence_type}{label};

    if(! exists $filter_info{'Evidence type'}{$ft_class_label}){
      $filter_info{'Evidence type'}{$ft_class_label} = 
        { count       => 0,
          description => $ftype_info->{$evidence_type}{long_name},
        };
    }
    $filter_info{'Evidence type'}{$ft_class_label}{count} += $count;
  }

  return \%filter_info;

  #Do we need to add an 'in_build' filter /data field?

}




#Need to bind param these as they come from URL parameters and are not tested

#Could move a lot of these to the BaseAdaptor if we have a valid_constraints config available
#and we use generic main_table approach for building constraint

#and reuse between adaptors if we use the _tables method to get the table syn
#This may mean contraints can be specified for classes which do not contain 
#the relevant fields.
#Allow this flexiblity or validate fields/constraint?
#Or implicit by location of contraint config, i.e. put it in the relevant 
#parent adaptors


#All these _constrain methods must return a valid constraint string, and a hashref of any other constraint config

sub _constrain_projects{
  my ($self, $egs) = @_;

  #enable query extension
  my $constraint_conf = {tables => [['input_set', 'inp'], ['experiment', 'e']]};
  
  
  if ( (ref($egs) ne 'ARRAY') ||
       scalar(@$egs) == 0 ) {
    throw('Must pass an arrayref of project names');
  }
  my @eg_ids;
  
  foreach my $eg (@$egs) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ExperimentalGroup', $eg);
    
    if (! $eg->is_project) {
      throw("You have passed an ExperimentalGroup which is not a project:\t".$eg->name);
    }
    
    push @eg_ids, $eg->dbID;
  }
  
  #Don't need to bind param this as we validate above
  my $constraint = ' fs.input_set_id=inp.input_set_id and inp.experiment_id=e.experiment_id AND '.
    'e.experimental_group_id IN ('.join(', ', @eg_ids).')';
  
  return ($constraint, $constraint_conf);
}


sub _constrain_evidence_types {
  my ($self, $etypes) = @_;

	my $constraint_conf = {tables => [['feature_type', 'ft']]};
     
  my %in_values = 
		(
		 'DNase1 & TFBS' => ['Open Chromatin', 
                         'Transcription Factor', 
                         'Transcription Factor Complex'],
     
		 'Hists & Pols'  => ['Histone',
                         'Polymerase'],
		);
      
  my @ft_classes;
      
  if ( (ref($etypes) ne 'ARRAY') ||
       scalar(@$etypes) == 0 ) {
		throw('Must pass an arrayref of evidence types');
  }
      
  foreach my $etype (@$etypes) {
    
		if (! exists $in_values{$etype}) {
		  throw("You have passed an invalid evidence type filter argument($etype)\n".
            "Please use one of the following:\t".join(' ,', keys(%in_values)));
		}
		push @ft_classes, @{$in_values{$etype}};
  }

  #Don't need to bind param this as we validate above
  my $constraint = ' fs.feature_type_id=ft.feature_type_id AND ft.class IN ("'.
		join('", "', @ft_classes).'")';
  
  return ($constraint, $constraint_conf);
}


# These following are duplicated in ResultSetAdaptor and potentially InputSetAdaptor
# Move to a new SetAdaptor? (not appropriate for DataSets/InputSubsets)

sub _constrain_cell_types {
  my ($self, $cts) = @_;

  my $constraint = ' fs.cell_type_id IN ('.
		join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $cts, 'dbID')}
        ).')';
     
  return ($constraint, {});  #{} = no futher contraint config
}


sub _constrain_feature_types {
  my ($self, $fts) = @_;
 

  my @tables = $self->_tables;
  my (undef, $syn) = @{$tables[0]};

  #Don't need to bind param this as we validate
  my $constraint = " ${syn}.feature_type_id IN (".
		join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fts, 'dbID')}).')';  
  
  return ($constraint, {});   #{} = not futher constraint conf
}

sub _constrain_analyses {
  my ($self, $anals) = @_;
  
  #Don't need to bind param this as we validate
  my $constraint = ' fs.analysis_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Analysis', $anals, 'dbID')}).')';
  
  return ($constraint, {});   #{} = not futher constraint conf
}

# add other fetch args
#type
#name
  





1;

__END__

#Methods to add?

#fetch_all_by_InputSet - would require input_set_id index
