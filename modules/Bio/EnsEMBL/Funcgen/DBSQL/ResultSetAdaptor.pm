#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ResultSetAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::ResultSetAdaptor - A database adaptor for fetching and
storing ResultSet objects.

=head1 SYNOPSIS

my $rset_adaptor = $db->get_ResultSetAdaptor();

my @rsets = @{$rset_adaptor->fetch_all_by_Experiment()};

=head1 DESCRIPTION

The ResultSetAdaptor is a database adaptor for storing and retrieving
ResultSet objects which encapsulate ResultFeatures defining individual points of 
experimental 'signal' data. This can be raw signal(Channel) or 
normalised(ExperimentalChip) data from an Array Experiment. A ResultSet can also 
encapsulate processed signal/read data(InputSet) from a sequencing Experiment e.g RPKM.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor;


use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::ResultSet;
use Bio::EnsEMBL::Funcgen::ResultFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(mean median);
use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor); #@ISA

# Table defs for use with dynamic query composition

use constant TRUE_TABLES => [[ 'result_set',        'rs' ],
                             [ 'result_set_input',  'rsi'],
                             [ 'dbfile_registry',   'dr' ]];

use constant TABLES => [[ 'result_set',        'rs' ],
                        [ 'result_set_input',  'rsi'],
                        [ 'dbfile_registry',   'dr' ]];


=head2 fetch_all_by_feature_class

  Arg [1]    : String - feature class i.e. 'result' or 'dna_methylation'.
  Arg [2]    : HASH of parameters (optional) containing contraint config e.g.

                   $result_set_adaptor->fetch_all_displayable_by_feature_class
                                           ('dna_methylation', 
                                             {'constraints' => 
                                               {
                                               cell_types     => [$cell_type],  #Bio::EnsEMBL::Funcgen::CellType
                                               #projects       => ['ENCODE'],
                                               feature_types  => [$ftype],      #Bio::EnsEMBL::Funcgen::FeatureType
                                               status         => 'DISPLAYABLE',
                                               } 
                                             });

  Example    : my @result_sets = @{$rs_adaptopr->fetch_all_by_feature_class('result')};
  Description: Retrieves ResultSet objects from the database based on result_set feature_class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::ResultSet objects
  Exceptions : Throws if type not defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_feature_class {
  my ($self, $fclass, $params) = @_;
  
  throw('Must provide a feature_set type') if(! defined $fclass);	
  my $sql = 'rs.feature_class = "'.$fclass.'"';
    
  #Deal with params constraints
  my $constraint = $self->compose_constraint_query($params);
  $sql .=  " AND $constraint " if $constraint;


  #Get result and reset true tables
  my $result = (defined $sql) ? $self->generic_fetch($sql) : [];
  $self->reset_true_tables;

  return $result;
}





=head2 fetch_all_linked_by_ResultSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::ResultSet
  Arg [2]    : Bio::EnsEMBL::Analysis
  Example    : my @rsets = @{$rset_adaptor->fetch_all_linked_by_ResultSet($rset)};
  Description: Retrieves a list of Bio::EnsEMBL::Funcgen::ResultSets which are linked
               to the supplied ResultSet (i.e. replicate relationships/import sets)
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ResultSet objects
  Exceptions : Throws if ResultSet not valid or stored
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_linked_by_ResultSet{
  my ($self, $rset) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);


  my $constraint = ' rsi.result_set_id in (SELECT distinct(result_set_id) from result_set_input where result_set_input_id in('.join(', ', @{$rset->result_set_input_ids}).')) ';
  
  my @tmp = @{$self->generic_fetch($constraint)};

  #Now remove query set
  my @linked_sets;

  map {push @linked_sets, $_ if $_->dbID != $rset->dbID} @tmp;

  return \@linked_sets;

}



=head2 fetch_all_by_Experiment_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Arg [2]    : Bio::EnsEMBL::Analysis
  Example    : my @rsets = @{$rset_adaptor->fetch_all_by_Experiment_Analysis($exp, $anal)};
  Description: Retrieves a list of Bio::EnsEMBL::Funcgen::ResultSets with the given Analysis from the Experiment
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ResultSet objects
  Exceptions : Throws if Analysis is not valid and stored
  Caller     : General
  Status     : At risk - to be merged into fetch_all_by_Experiment

=cut

sub fetch_all_by_Experiment_Analysis{
  my ($self, $exp, $analysis) = @_;

  if ( !($analysis && $analysis->isa("Bio::EnsEMBL::Analysis") && $analysis->dbID())) {
    throw("Need to pass a valid stored Bio::EnsEMBL::Analysis");
  }

  my $join = $self->_get_Experiment_join_clause($exp);
	
  return ($join) ?  $self->generic_fetch($join." AND rs.analysis_id=".$analysis->dbID) : [];
}

sub _get_Experiment_join_clause{
  my ($self, $exp) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exp);
  my $constraint;
  my @ecs = @{$exp->get_ExperimentalChips()};

  if (@ecs) { # We have an Array based experiment

    my $ec_ids = join(', ', (map $_->dbID, @ecs)); #get ' separated list of ecids
	
	
    my @chans = map @$_, (map $_->get_Channels(), @ecs);
    my $chan_ids = join(', ', (map $_->dbID(), @chans)); #get ' separated list of chanids
    #These give empty strings which are defined
    #This will not work for single IDs of 0, but this will never happen.
	
    if ($ec_ids && $chan_ids) {
      $constraint = '(((rsi.table_name="experimental_chip" AND rsi.table_id IN ('.$ec_ids.
        ')) OR (rsi.table_name="channel" AND rsi.table_id IN ('.$chan_ids.'))))';
      #This could probably be sped up using UNION
      #But result set is too small for cost of implementation
    } 
    elsif ($ec_ids) {
      $constraint = 'rsi.table_name="experimental_chip" AND rsi.table_id IN ('.$ec_ids.')';
    } 
    elsif ($chan_ids) {
      $constraint = 'rsi.table_name="channel" AND rsi.table_id IN ('.$chan_ids.')';
    }

  } 
  else {     #We have an InputSet Experiment
    #We could possibly have an experiment with an array and an input set
    #Currently nothing to stop this, but would most likely be loaded as separate experiments
    my $input_setids = join(', ', (map $_->dbID, @{$self->db->get_InputSetAdaptor->fetch_all_by_Experiment($exp)}));
    $constraint = 'rsi.table_name="input_set" AND rsi.table_id IN ('.$input_setids.')';
  }

  return $constraint;
}


=head2 fetch_all_by_Experiment

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my @rsets = @{$rset_adaptor->fetch_all_by_Experiment($exp)};
  Description: Retrieves a list of Bio::EnsEMBL::Funcgen::ResultSets from the Experiment
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ResultSet objects
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Experiment{
  my ($self, $exp, $analysis) = @_;

  #my $constraint = "ec.experiment_id=".$exp->dbID();
  #This was much easier with the more complicated default where join
  #should we reinstate and just have a duplication of cell/feature_types?
	
  my $join = $self->_get_Experiment_join_clause($exp);
  
  return ($join) ? $self->generic_fetch($join) : [];
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

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);	
  my $constraint = "rs.feature_type_id =".$ftype->dbID();
	
  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_name

  Arg [0]    : Mandatory string - ResultSet name
  Arg [1]    : Optional Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : Optional Bio::EnsEMBL::Funcgen::CellType
  Arg [3]    : Optional Bio::EnsEMBL::Analysis
  Example    : ($rset) = @{$rseta->fetch_all_by_name($exp->name().'_IMPORT')};
  Description: Retrieves ResultSets based on the name attribute
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::ResultSet objects
  Exceptions : Throws if no name provided or optional arguments are not valid   
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_name{
  my ($self, $name, $ftype, $ctype, $anal) = @_;

  if ( ! defined $name) {
    throw('Need to pass a ResultSet name');
  }
	
  my $constraint = 'rs.name = ?';
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);

  if ($ftype) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);
    $constraint .= ' AND rs.feature_type_id=?';
    $self->bind_param_generic_fetch($ftype->dbID, SQL_INTEGER);
  }

  if ($ctype) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType',    $ctype);
    $constraint .= ' AND rs.cell_type_id=?';
    $self->bind_param_generic_fetch($ctype->dbID, SQL_INTEGER);
  }

  if ($anal) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis',              $anal);
    $constraint .= ' AND rs.analysis_id=?';
    $self->bind_param_generic_fetch($anal->dbID, SQL_INTEGER);
  }

	
  return $self->generic_fetch($constraint);
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

  return ( @{$self->TABLES} );
}


=head2 _left_join

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the left join clasnames and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : Stable

=cut

#sub _left_join {
#  return (['dbfile_registry', '(rs.result_set_id=dr.table_id AND dr.table_name="result_set")'] , ['result_set_input','rs.result_set_id=rsi.result_set_id']);
#}
#Allows for absent result_set_input entries, in conjunction with omiting _default_where

sub _left_join {
  return (['dbfile_registry', '(rs.result_set_id=dr.table_id AND dr.table_name="result_set")']);
}


=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Stable

=cut

sub _columns {
	return qw(
            rs.result_set_id           rs.analysis_id
            rsi.table_name             rsi.result_set_input_id
            rsi.table_id               rs.name
            rs.cell_type_id            rs.feature_type_id
            rs.feature_class           dr.path
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
  Status     : Medium Risk

=cut

sub _default_where_clause {
  return 'rs.result_set_id = rsi.result_set_id';
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
  Status     : Medium Risk

=cut


sub _final_clause {
  #do not mess with this!
  return ' GROUP by rsi.result_set_input_id, rsi.result_set_id '.
    'ORDER BY rs.result_set_id, rs.cell_type_id, rs.feature_type_id';
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
  my ($sql, $table_name, $cc_id, $ftype_id, $ctype_id, $rf_set, $dbfile_dir,$feat_class);
  my $a_adaptor = $self->db->get_AnalysisAdaptor();
  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
  my $ct_adaptor = $self->db->get_CellTypeAdaptor(); 
  $sth->bind_columns(\$dbid, \$anal_id, \$table_name, \$cc_id, 
					 \$table_id, \$name, \$ctype_id, \$ftype_id, \$feat_class, \$dbfile_dir);
  
    #Need c/ftype cache here or rely on query cache?

  while ( $sth->fetch() ) {
    
    if(! $rset || ($rset->dbID() != $dbid)){
      
      push @rsets, $rset if $rset;

      $anal  = (defined $anal_id)  ? $a_adaptor->fetch_by_dbID($anal_id)   : undef;
      $ftype = (defined $ftype_id) ? $ft_adaptor->fetch_by_dbID($ftype_id) : undef;
      $ctype = (defined $ctype_id) ? $ct_adaptor->fetch_by_dbID($ctype_id) : undef;


      if(defined $dbfile_dir){
        $dbfile_dir = $self->build_dbfile_path($dbfile_dir);
      }
	 
      $rset = Bio::EnsEMBL::Funcgen::ResultSet->new
        (
         -DBID            => $dbid,
         -NAME            => $name,
         -ANALYSIS        => $anal,
         -TABLE_NAME      => $table_name,
         -FEATURE_TYPE    => $ftype,
         -CELL_TYPE       => $ctype,
         -FEATURE_CLASS   => $feat_class,
         -ADAPTOR         => $self,
         -DBFILE_DATA_DIR => $dbfile_dir,
        );
    }
    
  
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


#needs to move to a dbfile specific module

sub build_dbfile_path{
  my ($self, $subpath) = @_;

  my $full_path = $self->dbfile_data_root.'/'.$subpath;
  $full_path =~ s:/+:/:g;

  return $full_path;
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
  
  my $sth = $self->prepare('INSERT INTO result_set '.
                           '(analysis_id, name, cell_type_id, feature_type_id, feature_class) '.
                           'VALUES (?, ?, ?, ?, ?)');
  
  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();
  
 FEATURE: foreach my $rset (@rsets) {
    
    if ( ! (ref $rset && $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet') )) {
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
   

    my $ct_id = (defined $rset->cell_type)    ? $rset->cell_type->dbID    : undef;
    my $ft_id = (defined $rset->feature_type) ? $rset->feature_type->dbID : undef;

    $sth->bind_param(1, $rset->analysis->dbID,   SQL_INTEGER);
    $sth->bind_param(2, $rset->name,             SQL_VARCHAR);
    $sth->bind_param(3, $ct_id,                  SQL_INTEGER);
    $sth->bind_param(4, $ft_id,                  SQL_INTEGER);
    $sth->bind_param(5, $rset->feature_class,    SQL_VARCHAR);
    $sth->execute;
    
    $rset->dbID( $sth->{mysql_insertid} );
    $rset->adaptor($self);
    
    $self->store_states($rset);
    $self->store_chip_channels($rset);
    $self->store_dbfile_data_dir($rset) if $rset->dbfile_data_dir;
  }
  
  return \@rsets;
}





=head2 dbfile_data_root

  Arg[1]     : Optional String: Root path of dbfile data directory
  Example    : $rset_adaptor->dbfile_data_dir('/over-ride/path);
  Description: This allows the root path to be over-ridden. The default path
               is stored in the meta table as 'dbfile.data_root'. This can be
               over-ridden by the webcode by setting REGULATION_FILE_PATH in
               DEFAULTS.ini
  Returntype : String
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::DBAdaptor::ResultSet
  Status     : at risk - move this to SetAdaptor/FileAdaptor?

=cut

sub dbfile_data_root{
  my ($self, $root) = @_;

  if($root){
    $self->{dbfile_data_root} = $root;
  }
  elsif(! defined $self->{dbfile_data_root}){

    #Grab it from the meta table
    $self->{dbfile_data_root} = $self->db->get_MetaContainer->single_value_by_key('dbfile.data_root');
	
    #set to empty string if undef, so we don't keep querying the meta table
    $self->{dbfile_data_root} ||= '';
  }

  return $self->{dbfile_data_root};
}





=head2 store_dbfile_data_dir

  Arg[1]     : Bio::EnsEMBL::Funcgen::ResultSet
  Example    : $rset_adaptor->store_dbfile_data_dir;
  Description: Updater/Setter for the root dbfile data directory for this ResultSet
  Returntype : None
  Exceptions : Throws if ResultSet is not stored and valid
               Throws if dbfile_dat_dir not defined
  Caller     : Bio::EnsEMBL::Funcgen::Collector::ResultFeature
  Status     : at risk

=cut

#This is actually not always a dir
#leave for now until we use core API for this


sub store_dbfile_data_dir{
  my ($self, $rset) = @_;


  if ( ! $rset->is_stored($self->db()) ) {
    throw('ResultSet must be stored in the database before storing chip_channel entries');
  }

  my $data_dir = $rset->dbfile_data_dir;

  if(! $data_dir){
	throw('You need to pass a dbfile_data_dir argument to update the ResultSet');
  }
  
  #Check we have a record
  my $sql = 'SELECT path from dbfile_registry where table_name="result_set" and table_id=?';
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $data_dir, SQL_VARCHAR);
  $sth->execute;

  my ($db_dir) = $sth->fetch();

  if($db_dir &&
	 ($db_dir ne $data_dir)){#UPDATE
	$sql = 'UPDATE dbfile_registry set path=? where table_name="result_set" and table_id=?';
	$sth = $self->prepare($sql);
	$sth->bind_param(1, $data_dir,   SQL_VARCHAR);
	$sth->bind_param(2, $rset->dbID, SQL_INTEGER);
	$sth->execute;
  }
  elsif(! $db_dir){#STORE
	$sql = 'INSERT INTO  dbfile_registry(table_id, table_name, path) values(?, "result_set", ?)';
	$sth = $self->prepare($sql);
	$sth->bind_param(1, $rset->dbID, SQL_INTEGER);
	$sth->bind_param(2, $data_dir,   SQL_VARCHAR);
	
	$sth->execute;
  }
  
  return;
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
		INSERT INTO result_set_input (
			result_set_id, table_id, table_name
		) VALUES (?, ?, ?)
	");

  my $sth1 = $self->prepare("
		INSERT INTO result_set_input (
			result_set_input_id, result_set_id, table_id, table_name
		) VALUES (?, ?, ?, ?)
	");
  

  #Store and set all previously unstored table_ids
  foreach my $table_id(@{$rset->table_ids()}){
    my $cc_id = $rset->get_result_set_input_id($table_id);

    if(! defined $cc_id){
      $sth->bind_param(1, $rset->dbID(),       SQL_INTEGER);
      $sth->bind_param(2, $table_id,           SQL_INTEGER);
      $sth->bind_param(3, $rset->table_name(), SQL_VARCHAR);
      
      $sth->execute();

	  $cc_id = $sth->{'mysql_insertid'};
      $rset->add_table_id($table_id,  $sth->{'mysql_insertid'});
    }else{

	  #this should only store if not already stored for this rset
	  #this is because we may want to add chip_channels to a previously stored rset
	  my $sql = 'SELECT result_set_input_id from result_set_input where result_set_id='.$rset->dbID().
		" AND result_set_input_id=${cc_id}";
	  my ($loaded) = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

	  if(! $loaded){
		$sth1->bind_param(1, $cc_id,       SQL_INTEGER);
		$sth1->bind_param(2, $rset->dbID(),       SQL_INTEGER);
		$sth1->bind_param(3, $table_id,           SQL_INTEGER);
		$sth1->bind_param(4, $rset->table_name(), SQL_VARCHAR);
		$sth1->execute();#this could still fail is some one duplicates a result_set_id, table_id, table_name entry
	  }
	}
  }
  return $rset;
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
  Exceptions : Warns if not experimental_chip ResultSet
               Throws if no Slice passed
               Warns if 
  Caller     : general
  Status     : At risk

=cut


#Could we also have an optional net size for grouping ResultFeature into  Nbp pseudo ResultFeatures?


#This does not account for strandedness!!
###???
#Is this sensible?  Do we really want the full probe object alongside the ResultFeatures?
#strandedness?
#what about just the probe name?
#we could simplt add the probe_id to the ResultFeature
#This would prevent creating probe features for all of the features which do not have results for a given resultset
#This will mean the probe will still have to be individually created
#But we're only creating it for those which we require
#and we're now creating the lightweight ResultFeature instead of the ProbeFeature
#However, if we're dealing with >1 rset in a loop
#Then we'll be recreating the same ResultFeatures and probes for each set.
#We're already restricting the ProbeFeatures to those within the rset
#What we want is to get the score along side the ProbeFeature?
#But we want the probe name!!
#We really want something that will give Probe and ResultFeature
#Let's set the Probe as an optional ResultFeature attribute


sub fetch_ResultFeatures_by_Slice_ResultSet{
  my ($self, $slice, $rset, $ec_status, $with_probe) = @_;

  warn "Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor::fetch_ResultFeatures_by_Slice_ResultSet is now deprecated\n".
	"Please use the ResultFeatureAdaptor directly";

  return $self->db->get_ResultFeatureAdaptor->fetch_all_by_Slice_ResultSet($slice, $rset, $ec_status, $with_probe);

}



# Dynamic query contraint methods
# Most of these a re generic and redundant wrt FeatureSetAdaptor
# Could be moved to a SetAdaptor


sub _constrain_cell_types {
  my ($self, $cts) = @_;

  my @tables = $self->_tables;
  my (undef, $syn) = @{$tables[0]};

  my $constraint = " ${syn}.cell_type_id IN (".
		join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $cts, 'dbID')}
        ).')';
  
  #{} = no futher contraint config
  return ($constraint, {});
}


sub _constrain_feature_types {
  my ($self, $fts) = @_;
 
  my @tables = $self->_tables;
  my (undef, $syn) = @{$tables[0]};

  #Don't need to bind param this as we validate
  my $constraint = " ${syn}.feature_type_id IN (".
		join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fts, 'dbID')}).')';  
  
  #{} = not futher constraint conf
  return ($constraint, {});
}







### DEPRECATED METHODS


sub fetch_all_by_name_Analysis {#Deprecated in v69
  my ($self, $name, $anal) = @_;

  deprecate('To be removed in v71. Please use fetch_all_by_name');

  return $self->fetch_all_by_name($name, undef, undef, $anal);
}







1;

