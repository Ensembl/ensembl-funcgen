#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor - A database adaptor for fetching and
storing Funcgen Experiment objects.

=head1 SYNOPSIS

my $exp_a = $db->get_ExperimentAdaptor();
my $exp   = $exp_a->fetch_by_name($name);


=head1 DESCRIPTION

The ExperimentAdaptor is a database adaptor for storing and retrieving
Funcgen Experiment objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate warning );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;# For DBI :sql_types import;

use base qw( Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor );

=head2 fetch_by_name

  Arg [1]    : string - name of an Experiment
  Example    : my $exp = $exp_a->fetch_by_name('Exp-1');
  Description: Retrieves a named Experiment object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : Throws if no name defined or if more than one returned
  Caller     : General
  Status     : Medium risk

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw("Need to specify and experiment name argument") if  ! defined $name;

  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  my $result = $self->generic_fetch("e.name = ?");

  if (scalar @$result > 1) {
    throw("Experiment $name is not unique in the database, but only one result has been returned");
    #should have unique key of group_id and experiment_name
  }

  return $result->[0];
}


=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType object
  Example    : my @exps = @{$exp_adaptor->fetch_all_by_FeatureType($ftype)};
  Description: Retrieves Experiment objects from the database based on FeatureType
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_FeatureType {
  my $self   = shift;
  my $ftype  = shift;
  my $params = {constraints => {feature_types => [$ftype]}};
  return $self->generic_fetch($self->compose_constraint_query($params));
}



=head2 fetch_all_by_Epigenome

  Arg [1]    : Bio::EnsEMBL::Funcgen::Epigenome object
  Example    : my @exps = @{$exp_adaptor->fetch_all_by_Epigenome($epigenome)};
  Description: Retrieves Experiment objects from the database based on Epigenome
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_Epigenome {
  my $self       = shift;
  my $epigenome  = shift;
  my $params     = {constraints => {epigenomes => [$epigenome]}};
  return $self->generic_fetch($self->compose_constraint_query($params));
}


=head2 fetch_all_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Funcgen::Analysis
  Example    : my @exps = @{$exp_adaptor->fetch_all_by_Analysis($analysis)};
  Description: Retrieves Experiment objects from the database based on an Analysis
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Analysis {
  my ($self, $epigenome, $status) = @_;
  my $params = {constraints => {analyses => [$epigenome]}};
  $params->{constraints}{states} = [$status] if defined $status;
  my $results = $self->generic_fetch($self->compose_constraint_query($params));
  $self->reset_true_tables; #As we may have added status
  return $results;
}



=head2 get_all_experiment_names

  Arg [1]    : (optional) boolean - flag to denote whether experiment is flagged for web display
  Example    : my @names = @{$exp_a->get_all_experiment_names()};
  Description: Retrieves names of all experiments.
  Returntype : ARRAYREF
  Exceptions : none
  Caller     : General
  Status     : At Risk - rename fetch?

=cut

sub get_all_experiment_names{
  my $self        = shift;
  my $displayable = shift;
  my $sql         = 'SELECT e.name FROM experiment e';

  if($displayable){
    $sql .= ', status s, status_name sn WHERE e.experiment_id = s.table_id AND '.
      's.table_name="experiment" AND s.status_name_id = sn.status_name_id'.
      ' and sn.name="DISPLAYABLE"';
  }

  return $self->db->dbc->db_handle->selectcol_arrayref($sql);
}


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
  return (['experiment', 'e']);
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
    return qw( e.experiment_id e.name e.experimental_group_id
        e.control_id e.is_control
        e.feature_type_id e.epigenome_id e.archive_id );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
    @result,       $exp_id,     $name,
    $group_id,     $control_id, $is_control,
    $epigenome_id, $ft_id,      $archive_id,
  );

  my $eg_adaptor   = $self->db->get_ExperimentalGroupAdaptor;
  my $epi_adaptor  = $self->db->get_EpigenomeAdaptor;
  my $ft_adaptor   = $self->db->get_FeatureTypeAdaptor;
  my $exp_adaptor  = $self->db->get_ExperimentAdaptor;

  $sth->bind_columns( 
    \$exp_id, \$name, \$group_id, \$control_id,
    \$is_control, \$ft_id, \$epigenome_id, \$archive_id
  );

  my (%ftypes, %epigenomes, $control);

  while ( $sth->fetch() ) {

    my $group = $eg_adaptor->fetch_by_dbID($group_id);#cache these in ExperimentalGroupAdaptor
    my $current_epigenome = undef;

    if(defined $epigenome_id && ! exists $epigenomes{$epigenome_id}) {
      $epigenomes{$epigenome_id} = $epi_adaptor->fetch_by_dbID($epigenome_id);
      if(! defined $epigenomes{$epigenome_id}){
        throw("Could not fetch linked Epigenome (dbID: $epigenome_id) for Experiment:\t$name");
      }
      $current_epigenome = $epigenomes{$epigenome_id};
    }

    if(! exists $ftypes{$ft_id}) {
      $ftypes{$ft_id} = $ft_adaptor->fetch_by_dbID($ft_id);
      if(! defined $ftypes{$ft_id}) {
        throw("Could not fetch linked FeatureType (dbID: $ft_id) for Experiment:\t$name");
      }
    }

    if($control_id) {
      $control = $exp_adaptor->fetch_by_dbID($control_id);
      if(! $control) {
        throw("Could not fetch linked control Experiment (dbID: $control_id) for Experiment:\t$name");
      }
    }

    push @result, Bio::EnsEMBL::Funcgen::Experiment->new(
      -DBID                => $exp_id,
      -ADAPTOR             => $self,
      -NAME                => $name,
      -FEATURE_TYPE        => $ftypes{$ft_id},
      -EPIGENOME           => $current_epigenome,
      -ARCHIVE_ID          => $archive_id,
      -EXPERIMENTAL_GROUP  => $group,
      -CONTROL             => $control,
      -IS_CONTROL          => $is_control,
    );
  }
  return \@result;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::Experiment objects
  Example    : $oaa->store($exp1, $exp2, $exp3);
  Description: Stores given Experiment objects in the database.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : Throws is group not present in DB
               Throws if object is not a Bio::EnsEMBL::Funcgen::Experiment
               Throws if object is already present in the DB but has no dbID
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @exps = @_;

    my $sth = $self->prepare(
        'INSERT INTO experiment(name, experimental_group_id, control_id,
                            is_control,feature_type_id, epigenome_id, archive_id)
                            VALUES (?, ?, ?, ?, ?, ?, ?)'
    );




  foreach my $exp (@exps) {
    assert_ref($exp, 'Bio::EnsEMBL::Funcgen::Experiment');
  
	  if (! $exp->is_stored($self->db)){
      #Test object attrs are stored
      my $exp_group = $exp->experimental_group;
      $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ExperimentalGroup', $exp_group);
      $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType',       $exp->feature_type);
      $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Epigenome',         $exp->epigenome);
    if ( $exp->get_control ) {
        $self->db->is_stored_and_valid( 'Bio::EnsEMBL::Funcgen::Experiment',
            $exp->get_control );
    }


		  #Validate doesn't exist aleady
		
		  if($self->fetch_by_name($exp->name)){
  		  throw('Experiment ['.$exp->name.'] already exists in the database.'.
			   "\nTo reuse/update this Experiment you must retrieve it using the ExperimentAdaptor");
		  }
  
  		$sth->bind_param(1,  $exp->name,                      SQL_VARCHAR);
      $sth->bind_param(2,  $exp_group->dbID,                SQL_INTEGER);

      my $control_id = undef;
      if ($exp->get_control){
        $control_id = $exp->get_control->dbID;
      }
      
      $sth->bind_param(3,  $control_id,                     SQL_INTEGER);
  		$sth->bind_param(4,  $exp->is_control,                SQL_TINYINT);
      $sth->bind_param(5,  $exp->feature_type->dbID,        SQL_INTEGER); 
      $sth->bind_param(6,  $exp->epigenome->dbID,           SQL_INTEGER);
      $sth->bind_param(7, $exp->archive_id,                SQL_VARCHAR); 
  		$sth->execute();
  		$exp->dbID($self->last_insert_id);
  		$exp->adaptor($self);
	  }
	  else{ #assume we want to update the states
		  warning('If you are trying to update the states, you may want to use $exp->adaptor->store_states($exp)');
		  $self->store_states($exp);
	  }
	}

    return \@exps;
}



=head2 fetch_source_label_by_experiment_id

  Args       : Int - experiment_id
  Example    : my $source_label = $fset->source_label;
  Description: Retrieves the source label this FeatureSet, used in zmenus
  Returntype : String - Space separated if more than 1 label.
  Exceptions : 
  Caller     : FeatureSet/ResultSet::source_label, ultimately the webcode
  Status     : at risk

=cut

#This shortcuts having to create the Experiment object
#which probably saves 3 trips to the DB and some redundant processing

sub fetch_source_label_by_experiment_id{
  my $self   = shift;
  my $exp_id = shift or throw('Must provide an experiment_id argument');

  my $sql = 'SELECT e.archive_id,  eg.name, eg.is_project from experiment e '.
    'LEFT JOIN experimental_group eg using(experimental_group_id) where e.experiment_id=?';

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $exp_id, SQL_INTEGER);
 
  if(! eval{ $sth->execute; 1}){
    throw("Failed to fetch_source_label_by_experiment_id, SQL:\n$sql\n$@");
  }
  
  my ($archive_id,  $eg_name, $is_project) = $sth->fetchrow_array;
  $sth->finish;

  # Handle multiple SRX IDs, just in case the submitters
  # archive the reps as separate experiments :(
  my @source_labels;
  @source_labels = split/,/, $archive_id if defined $archive_id;
  push @source_labels, $eg_name if $is_project;  #Append project name
  
  return join(q{ }, # Single space
              @source_labels);
}




=head2 fetch_mage_xml_by_Experiment

  Args       : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my $xml = $exp_adaptor->fetch_mage_xml_by_Experiment($exp);
  Description: Gets the MAGE XML for this experiment
  Returntype : string
  Exceptions : throws if arg is not a valid stored Experiment
  Caller     : general
  Status     : at Risk

=cut

sub fetch_mage_xml_by_Experiment{
  my ($self, $exp) = @_;

  if(!($exp and $exp->isa('Bio::EnsEMBL::Funcgen::Experiment') && $exp->dbID())){
	throw('You must provide a valid stored Bio::EnsEMBL::Funcgen::Experiment');
  }

  return if ! $exp->mage_xml_id();

  my $sql = 'SELECT xml FROM mage_xml WHERE mage_xml_id='.$exp->mage_xml_id;

  return $self->db->dbc->db_handle->selectall_arrayref($sql)->[0];
}

=head2 fetch_mage_xml_by_experiment_name

  Args       :
  Example    : my $xml = $exp_adaptor->fetch_mage_xml_by_Experiment($exp);
  Description: Gets the MAGE XML for this experiment
  Returntype : string
  Exceptions : throws if no arg passed
  Caller     : general
  Status     : at Risk

=cut

sub fetch_mage_xml_by_experiment_name{
  my ($self, $exp_name) = @_;

  if(! defined $exp_name){
	throw('You must provide an Experiment name argument');
  }

  my $sql = 'SELECT mx.xml FROM mage_xml mx, experiment e WHERE e.name="'.$exp_name.'" and e.mage_xml_id=mx.mage_xml_id';

  return $self->db->dbc->db_handle->selectall_arrayref($sql)->[0];
}



=head2 update_mage_xml_by_Experiment

  Args       : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my $xml = $exp_adaptor->fetch_mage_xml_by_Experiment($exp);
  Description: Gets the MAGE XML for this experiment
  Returntype : string
  Exceptions : throws if arg is not a valid stored Experiment
  Caller     : general
  Status     : at Risk

=cut

sub update_mage_xml_by_Experiment{
  my ($self, $exp) = @_;


  if(!($exp and $exp->isa('Bio::EnsEMBL::Funcgen::Experiment'))){
	throw('You must provide a valid Bio::EnsEMBL::Funcgen::Experiment');
  }

  if($exp->mage_xml_id()){
	#potentially calling dbID on a un-stored obj, implicit that it
	warn('Overwriting mage_xml entry for Experiment:    '.$exp->name);
	my $sql = "UPDATE mage_xml set xml='".$exp->mage_xml()."'";
	$self->db->dbc->do($sql);

  }else{
	my $sql = "INSERT INTO mage_xml (xml) VALUES('".$exp->mage_xml()."')";
	#need to get a statement handle to retrieve insert id
	my $sth = $self->prepare($sql);
	$sth->execute();
	$exp->mage_xml_id($self->last_insert_id);

	$sql = "UPDATE experiment set mage_xml_id=".$exp->mage_xml_id()." where experiment_id =".$exp->dbID();
	$sth = $self->prepare($sql);
	$sth->execute();
  }

  return $exp;
}


## GENERIC CONSTRAIN METHODS ###

#All these _constrain methods must return a valid constraint string, and a hashref of any other constraint config

#Need to bind param any of these which come from URL parameters and are not tested

sub _constrain_epigenomes {
  my ($self, $epigenomes) = @_;

  my $constraint = $self->_table_syn.'.epigenome_id IN ('.
        join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::Epigenome', $epigenomes, 'dbID')}
        ).')';

  #{} = no futher contraint config
  return ($constraint, {});
}


sub _constrain_feature_types {
  my ($self, $fts) = @_;

  #Don't need to bind param this as we validate
  my $constraint = $self->_table_syn.'.feature_type_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fts, 'dbID')}).')';

  #{} = no futher constraint conf
  return ($constraint, {});
}


sub _constrain_analyses {
  my ($self, $anals) = @_;

  #Don't need to bind param this as we validate
  my $constraint = $self->_table_syn.'.analysis_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Analysis', $anals, 'dbID')}).')';

  return ($constraint, {});   #{} = no futher constraint conf
}


#This will not support concat'd archive IDs

sub _constrain_archive_ids {
  my ($self, $archive_ids) = @_;

  if ( (ref($archive_ids) ne 'ARRAY') || scalar(@$archive_ids) == 0 ) {
    throw('Must pass an arrayref of archive IDs');
  }
  my $constraint;

  #{} = not futher constraint conf
  return (' e.archive_id IN '.join(', ', (map {uc($_)} @$archive_ids)) , {});
}


1;

