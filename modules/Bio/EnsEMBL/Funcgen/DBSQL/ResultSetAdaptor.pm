#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ResultSetAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::ResultSetAdaptor - A database adaptor for fetching and
storing ResultSet objects.

=head1 SYNOPSIS

my $rset_adaptor = $db->get_ResultSetAdaptor();

my @rsets = @{$rset_adaptor->fetch_all_by_Experiment};

=head1 DESCRIPTION

The ResultSetAdaptor is a database adaptor for storing and retrieving
ResultSet objects. A ResultSet represents low level analysis of a discrete set 
of data i.e. alignments of a particular technical replicate or a merged set of
replicates. The underlying data are stored in flat files  which can be accessed 
via the relevant library or command line tool. 

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception         qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::ResultSet;
use Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor; #DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor);


=head2 fetch_all_by_feature_class

  Arg [1]    : String - feature class i.e. 'result' or 'dna_methylation'.
  Arg [2]    : HASH of parameters (optional) containing contraint config e.g.

                   $result_set_adaptor->fetch_all_displayable_by_feature_class
                                           ('dna_methylation',
                                             {'constraints' =>
                                               {
                                               epigenomes     => [$epigenome],  #Bio::EnsEMBL::Funcgen::Epigenome
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

  my $constraint = 'rsi.result_set_id IN ( '.
    'SELECT DISTINCT(result_set_id) FROM result_set_input '.
    'WHERE result_set_input_id IN('.join(', ', @{$rset->result_set_input_ids}).')) ';
  my @tmp = @{$self->generic_fetch($constraint)};

  # Now remove query set
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
#   my @ecs = @{$exp->get_ExperimentalChips()};
# 
#   if (@ecs) { # We have an Array based experiment
# 
#     my $ec_ids = join(', ', map {$_->dbID} @ecs); #get ' separated list of ecids
# 
# 
#     my @chans = map { @$_ } (map {$_->get_Channels} @ecs);
#     my $chan_ids = join(', ', map {$_->dbID} @chans); #get ' separated list of chanids
#     #These give empty strings which are defined
#     #This will not work for single IDs of 0, but this will never happen.
# 
#     if ($ec_ids && $chan_ids) {
#       $constraint = '(((rsi.table_name="experimental_chip" AND rsi.table_id IN ('.$ec_ids.
#         ')) OR (rsi.table_name="channel" AND rsi.table_id IN ('.$chan_ids.'))))';
#       #This could probably be sped up using UNION
#       #But result set is too small for cost of implementation
#     }
#     elsif ($ec_ids) {
#       $constraint = 'rsi.table_name="experimental_chip" AND rsi.table_id IN ('.$ec_ids.')';
#     }
#     elsif ($chan_ids) {
#       $constraint = 'rsi.table_name="channel" AND rsi.table_id IN ('.$chan_ids.')';
#     }
# 
#   }
#   else {     #We have an InputSet/InputSubset Experiment
    my $setids = join(', ', map {$_->dbID} 
      @{$self->db->get_InputSubsetAdaptor->fetch_all_by_Experiments([$exp])});
    $constraint = "rsi.table_name='input_subset' AND rsi.table_id IN (${setids})";
#   }

  return $constraint;
}

#todo deprecate above method in favour of fetch_all_by_Experiment

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
  my $join = $self->_get_Experiment_join_clause($exp);
  return ($join) ? $self->generic_fetch($join) : [];
}

=head2 fetch_all_by_name

  Arg [0]    : Mandatory string - ResultSet name
  Arg [1]    : Optional Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : Optional Bio::EnsEMBL::Funcgen::Epigenome
  Arg [3]    : Optional Bio::EnsEMBL::Analysis
  Example    : ($rset) = @{$rseta->fetch_all_by_name($exp->name().'_IMPORT')};
  Description: Retrieves ResultSets based on the name attribute
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::ResultSet objects
  Exceptions : Throws if no name provided or optional arguments are not valid
  Caller     : General
  Status     : At Risk

=cut

#This not longer handles the unique key as we don't have feature_class here
#this is to move to dbfile_registry.format

sub fetch_all_by_name{
  my ($self, $name, $ftype, $epigenome, $anal) = @_;

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

  if ($epigenome) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Epigenome',    $epigenome);
    $constraint .= ' AND rs.epigenome_id=?';
    $self->bind_param_generic_fetch($epigenome->dbID, SQL_INTEGER);
  }

  if ($anal) {
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis',              $anal);
    $constraint .= ' AND rs.analysis_id=?';
    $self->bind_param_generic_fetch($anal->dbID, SQL_INTEGER);
  }


  return $self->generic_fetch($constraint);
}

# sub fetch_by_Epigenome_FeatureType {
# 
#   my ($self, $epigenome, $feature_type) = @_;
#   
#   my @constraint;
#   
#   if ($feature_type) {
#     $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $feature_type);
#     push @constraint, 'rs.feature_type_id=?';
#     $self->bind_param_generic_fetch($feature_type->dbID, SQL_INTEGER);
#   }
#   if ($epigenome) {
#     $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Epigenome',    $epigenome);
#     push @constraint, 'rs.epigenome_id=?';
#     $self->bind_param_generic_fetch($epigenome->dbID, SQL_INTEGER);
#   }
#   
#   my $constraint = join ' and ', @constraint;
#   
#   return $self->generic_fetch($constraint);
# }

=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  my $self = shift;

  return ([ 'result_set',        'rs' ],
          [ 'result_set_input',  'rsi'],
          [ 'dbfile_registry',   'dr' ]);
#           [ 'dbfile_registry',   'dr_bigwig' ],
#           [ 'dbfile_registry',   'dr_bam'    ]);
          ;
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
#   return ([ 'result_set_input', '(rs.result_set_id=rsi.result_set_id)' ],
#           ['dbfile_registry', '(rs.result_set_id=dr.table_id AND dr.table_name="result_set" and dr.file_type="BIGWIG")']);

  my $self = shift;
  
  my $file_type = $self->{file_type};
  
  $file_type = "BIGWIG" unless($self->{file_type});

  return ([ 'result_set_input', '(rs.result_set_id=rsi.result_set_id)' ],
          ['dbfile_registry', '(rs.result_set_id=dr.table_id AND dr.table_name="result_set" and dr.file_type="' . $file_type . '")']);
#   return ([ 'result_set_input', '(rs.result_set_id=rsi.result_set_id)' ],
#           ['dr_bigwig', '(rs.result_set_id=dr_bigwig.table_id AND dr_bigwig.table_name="result_set" and dr_bigwig.file_type="BIGWIG")'],
#           ['dr_bam', '(rs.result_set_id=dr_bam.table_id AND dr_bam.table_name="result_set" and dr_bam.file_type="BAM")']);
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
            rs.epigenome_id            rs.feature_type_id
            rs.feature_class           dr.path
            rs.experiment_id
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
  return '';
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
    'ORDER BY rs.result_set_id, rs.epigenome_id, rs.feature_type_id';
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
  my $self = shift;
  my $sth  = shift;
  my (@rsets, $rset, $dbid, $anal_id, $anal, $ftype, $epigenome, $table_id, $experiment);
  my ($sql, $table_name, $cc_id, $ftype_id, $epigenome_id, $dbfile_path);
  my ($name, $rep, $feat_class);
  my $experiment_id;
  my $a_adaptor  = $self->db->get_AnalysisAdaptor;
  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor;
  my $epi_adaptor = $self->db->get_EpigenomeAdaptor;
  my $experiment_adaptor = $self->db->get_ExperimentAdaptor;
  $sth->bind_columns(\$dbid, \$anal_id, \$table_name, \$cc_id, \$table_id,
                     \$name, \$epigenome_id, \$ftype_id, \$feat_class, \$dbfile_path, \$experiment_id);

  while ( $sth->fetch ) {

    if( (! defined $rset) || ($rset->dbID != $dbid) ){
      push @rsets, $rset if $rset;
      $anal       = (defined $anal_id)       ? $a_adaptor->fetch_by_dbID($anal_id)   : undef;
      $ftype      = (defined $ftype_id)      ? $ft_adaptor->fetch_by_dbID($ftype_id) : undef;
      $epigenome  = (defined $epigenome_id)  ? $epi_adaptor->fetch_by_dbID($epigenome_id) : undef;
      $experiment = (defined $experiment_id) ? $experiment_adaptor->fetch_by_dbID($experiment_id) : undef;
    
#       if(defined $dbfile_path){
#         ($dbfile_path = $self->dbfile_data_root.'/'.$dbfile_path) =~ s:/+:/:g;
#       }

      $rset = Bio::EnsEMBL::Funcgen::ResultSet->new
        (-DBID            => $dbid,
         -NAME            => $name,
         -ANALYSIS        => $anal,
         -TABLE_NAME      => $table_name,
         -FEATURE_TYPE    => $ftype,
         -EPIGENOME       => $epigenome,
         -FEATURE_CLASS   => $feat_class,
         -ADAPTOR         => $self,
         -DBFILE_PATH     => $dbfile_path,
         -EXPERIMENT      => $experiment,
         );
    }

    $rset->_add_table_id($table_id, $cc_id);
    
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

sub store {
  my ($self, @rsets) = @_;
  scalar(@rsets) || throw("Must provide a list of ResultSet objects");

  my $sth = $self->prepare('INSERT INTO result_set '.
                           '(analysis_id, name, epigenome_id, feature_type_id, feature_class) '.
                           'VALUES (?, ?, ?, ?, ?)');
  my $db = $self->db;

 FEATURE: foreach my $rset (@rsets) {

    if ( ! (ref $rset && $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet') )) {
      throw('Must be an ResultSet object to store. Passed: ' . ref($rset) );
    }

    if ( $rset->is_stored($db) ) {
      throw('ResultSet [' . $rset->dbID . '] is already stored in the database'.
        "\nResultSetAdaptor does not yet accomodate updating ResultSets");
    }

    $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $rset->analysis);

    my $epi_id = (defined $rset->epigenome)    ? $rset->epigenome->dbID    : undef;
    my $ft_id = (defined $rset->feature_type) ? $rset->feature_type->dbID : undef;

    $sth->bind_param(1, $rset->analysis->dbID,   SQL_INTEGER);
    $sth->bind_param(2, $rset->name,             SQL_VARCHAR);
    $sth->bind_param(3, $epi_id,                 SQL_INTEGER);
    $sth->bind_param(4, $ft_id,                  SQL_INTEGER);
    $sth->bind_param(5, $rset->feature_class,    SQL_VARCHAR);

    if(! eval {$sth->execute; 1}){
      throw("Failed to store $rset ".$rset->name."\n$@");
    }

    $rset->dbID( $self->last_insert_id );
    $rset->adaptor($self);
#     $self->store_states($rset);
    $self->store_chip_channels($rset);
    $self->store_dbfile_path($rset, 'DIR') if $rset->dbfile_path;
  }

  return \@rsets;
}

=head2 store_dbfile_path

  Arg[1]     : Bio::EnsEMBL::Funcgen::ResultSet
  Arg[2]     : File type as a string
  Example    : $rset_adaptor->store_dbfile_path($result_set, 'BAM');
  Description: Updater/Setter for the root dbfile data directory for this ResultSet
  Returntype : None
  Exceptions : Throws if ResultSet is not stored and valid
               Throws if dbfile_data_dir not defined
               Throws if unable to store dbfile_data_dir
  Caller     : Bio::EnsEMBL::Funcgen::Hive::BigWigWriter
  Status     : at risk

=cut

sub store_dbfile_path {
  my $self = shift;
  my $rset = shift;
  my $file_type = shift;
  
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);

  my $path = $rset->dbfile_path;
  if (! defined $path) {
    throw('ResultSet::dbfile_path attribute is not set') 
  }
  if (! defined $file_type) {
    throw('file_type parameter has not been set!');
  }
  my $root = $self->dbfile_data_root;
  if (! $root) {
    throw('It is unsafe to store_dbfile_path without setting the dbfile_data_root first');
  }
  
  use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_checksum );
  my $md5sum = generate_checksum($rset->dbfile_path);
  
  $path =~ s/$root//;

  #Check we have a record
  my $rset_id  = $rset->dbID;
  my $db_path  = $self->_fetch_dbfile_path($rset_id, $file_type);

  if($db_path &&
	 ($db_path ne $path)) {  # UPDATE
    # Really should have rolled this back prior to this point
    my $sql = 'UPDATE dbfile_registry set path=? where table_name="result_set" and table_id=? and file_type=? and md5sum=?';
    my $sth = $self->prepare($sql);
    $sth->bind_param(1, $path,    SQL_VARCHAR);
    $sth->bind_param(2, $rset_id, SQL_INTEGER);
    $sth->bind_param(3, $file_type);
    $sth->bind_param(4, $md5sum,  SQL_VARCHAR);

    if(! eval {$sth->execute; 1}) {
      throw('Failed to update dbfile_data_dir for '.$rset->name."\n$@");
    }
  }
  elsif(! defined $db_path) {  # STORE

    my $sql = 'INSERT INTO dbfile_registry(table_id, table_name, path, file_type, md5sum) values(?, "result_set", ?, ?, ?)';
    my $sth = $self->prepare($sql);
    $sth->bind_param(1, $rset_id, SQL_INTEGER);
    $sth->bind_param(2, $path,    SQL_VARCHAR);
    $sth->bind_param(3, $file_type);
    $sth->bind_param(4, $md5sum,  SQL_VARCHAR);

    if(! eval {$sth->execute; 1}) {
      my $err = $@;
      #This could be a race condition if we have parallel writes going on
      #Attempt to validate stored value is same, else fail
      $db_path = $self->_fetch_dbfile_path($rset_id, $file_type);

      if(defined $db_path) {
        if($db_path ne $path) {
          throw('Failed to store dbfile_data_dir table '.$rset->name.
            "\n'Racing' process stored a differing value:\n\t$path\n\tvs\n\t$db_path\n$err");
        }  # else this was a race condition
      }
      else{
        throw('Failed to store dbfile_path for '.$rset->name."\n$err");
      }
    }
  }
  return;
}

=head2 dbfile_data_root
=cut
sub dbfile_data_root {
  my $self = shift; 
  return $self->db->dbfile_data_root(@_);
}

#This is only used for validation in store_dbfile_path

sub _fetch_dbfile_path {
  my $self    = shift;
  my $rset_id = shift;
  my $file_type = shift;
  
  my $sql  = qq(SELECT path from dbfile_registry where table_name="result_set" and table_id=$rset_id and file_type="$file_type");
  
  my $sth  = $self->prepare($sql);

  if(! eval {$sth->execute; 1} ){
    throw("Failed to fetch dbfile_registry using:\n$sql (dbID=$rset_id)\n$@");
  }

  my $path;

  if(my $row_ref = $sth->fetch){
    $path = $row_ref->[0];
  }

  return $path;
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

  my $sth = $self->prepare('INSERT INTO result_set_input (result_set_id, table_id, table_name)'.
    ' VALUES (?, ?, ?)');

  my $sth1 = $self->prepare('INSERT INTO result_set_input '.
    '(result_set_input_id, result_set_id, table_id, table_name) VALUES (?, ?, ?, ?)');

  #Store and set all previously unstored table_ids
  foreach my $table_id(@{$rset->table_ids()}){
    my $cc_id = $rset->get_result_set_input_id($table_id);

    if(! defined $cc_id){
      $sth->bind_param(1, $rset->dbID,       SQL_INTEGER);
      $sth->bind_param(2, $table_id,         SQL_INTEGER);
      $sth->bind_param(3, $rset->table_name, SQL_VARCHAR);
      $sth->execute;
      $rset->_add_table_id($table_id, $self->last_insert_id);
    }
    else{
      #this should only store if not already stored for this rset
      #this is because we may want to add chip_channels to a previously stored rset
      my $sql = 'SELECT result_set_input_id from result_set_input where result_set_id='.$rset->dbID.
       " AND result_set_input_id=${cc_id}";
      my $loaded = $self->db->dbc->db_handle->selectcol_arrayref($sql);

      if(scalar @$loaded){
        $sth1->bind_param(1, $cc_id,            SQL_INTEGER);
        $sth1->bind_param(2, $rset->dbID,       SQL_INTEGER);
        $sth1->bind_param(3, $table_id,         SQL_INTEGER);
        $sth1->bind_param(4, $rset->table_name, SQL_VARCHAR);
        $sth1->execute();#this could still fail is some one duplicates a result_set_id, table_id, table_name entry
      }
    }
  }
  return $rset;
}


### GENERIC CONSTRAIN METHODS ###
# See Base/SetAdaptor

1;

