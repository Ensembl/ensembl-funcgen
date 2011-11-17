#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor - A database adaptor for fetching and
storing Funcgen Experiment objects.

=head1 SYNOPSIS

my $exp_a = $db->get_ExperimentAdaptor();
my $exp   = $exp_a->fetch_by_name($name);


=head1 DESCRIPTION

The ExperimentAdaptor is a database adaptor for storing and retrieving
Funcgen Experiment objects.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


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
  
  throw("Need to specify and experiment name argument") if (! defined $name);
  
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  my $result = $self->generic_fetch("e.name = ?");
  
  if (scalar @$result > 1) {
    throw("Experiment $name is not unique in the database, but only one result has been returned");
    #should have unique key of group_id and experiment_name
  } 
  return $result->[0];
}


=head2 fetch_by_archive_id

  Arg [1]    : String - Archive ID to a public repository (ENA) e.g. SRX00381237
  Example    : my $exp = $exp_a->fetch_by_archive_id('SRX00381237');
  Description: Retrieves an Experiment via it's archive ID.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : Throws if archive ID not defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_archive_id {
  my $self = shift;
  my $id = shift;
  
  throw("Need to specify and archive id argument") if (! defined $id);

  $self->bind_param_generic_fetch($id, SQL_VARCHAR);
  my $result = $self->generic_fetch("e.archive_id = ?");
 
  return $result->[0];
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
  my ($self, $displayable) = @_;


  my ($constraint);

  my $sql = "SELECT e.name FROM experiment e";
  $sql .= ", status s WHERE e.experiment_id =\"s.table_id\" AND s.table_name=\"experiment\" AND s.state=\"DISPLAYABLE\"" if($displayable);
  
  return $self->db->dbc->db_handle->selectcol_arrayref($sql);
}

=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At risk

=cut

sub _tables {
	my $self = shift;
	
	return ['experiment', 'e'];
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
	
	return qw(e.experiment_id e.name e.experimental_group_id e.date e.primary_design_type e.description e.archive_id e.data_url e.mage_xml_id);
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
	
	my (@result, $exp_id, $name, $group_id, $p_design_type, $date, $description, $archive_id, $data_url, $xml_id);
	
	my $eg_adaptor = $self->db->get_ExperimentalGroupAdaptor();

	$sth->bind_columns(\$exp_id, \$name, \$group_id, \$date, \$p_design_type, \$description, \$archive_id, \$data_url, \$xml_id);
	

	while ( $sth->fetch() ) {

	  my $group = $eg_adaptor->fetch_by_dbID($group_id);#cache these in ExperimentalGroupAdaptor

	  my $exp = Bio::EnsEMBL::Funcgen::Experiment->new(
							   -DBID                => $exp_id,
							   -ADAPTOR             => $self,
							   -NAME                => $name,
							   -EXPERIMENTAL_GROUP  => $group,
							   -DATE                => $date,
							   -PRIMARY_DESIGN_TYPE => $p_design_type,
							   -DESCRIPTION         => $description,
							   -ARCHIVE_ID        => $archive_id,
							   -DATA_URL            => $data_url,
							   -MAGE_XML_ID         => $xml_id,
													  );
	  
	  push @result, $exp;
	  
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
    my @args = @_;
	
	my ($s_exp);
   	
	my $sth = $self->prepare('INSERT INTO experiment
                                 (name, experimental_group_id, date, primary_design_type, description, archive_id, data_url, mage_xml_id)
                                 VALUES (?, ?, ?, ?, ?, ?, ?, ?)');

    foreach my $exp (@args) {
	  throw('Can only store Experiment objects') 	if ( ! $exp->isa('Bio::EnsEMBL::Funcgen::Experiment'));
	  
	  if (!( $exp->dbID() && $exp->adaptor() == $self )){
		
		
		my $exp_group = $exp->experimental_group;
		$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ExperimentalGroup', $exp_group);

		
		$s_exp = $self->fetch_by_name($exp->name());#validate on group too!
		throw("Experimental already exists in the database with dbID:".$s_exp->dbID().
			  "\nTo reuse/update this Experimental you must retrieve it using the ExperimentalAdaptor".
			  "\nMaybe you want to use the -recover option?") if $s_exp;
		
		
		$exp = $self->update_mage_xml_by_Experiment($exp) if(defined $exp->mage_xml());
			
		$sth->bind_param(1, $exp->name(),                     SQL_VARCHAR);
		$sth->bind_param(2, $exp->experimental_group()->dbID, SQL_INTEGER);
		$sth->bind_param(3, $exp->date(),                     SQL_VARCHAR);#date?
		$sth->bind_param(4, $exp->primary_design_type(),      SQL_VARCHAR);
		$sth->bind_param(5, $exp->description(),              SQL_VARCHAR);
		$sth->bind_param(6, $exp->archive_id(),               SQL_VARCHAR);
		$sth->bind_param(7, $exp->data_url(),                 SQL_VARCHAR);
		$sth->bind_param(8, $exp->mage_xml_id(),              SQL_INTEGER);
		
		$sth->execute();
		$exp->dbID($sth->{'mysql_insertid'});
		$exp->adaptor($self);
		
	  }
	  else{
		#assume we want to update the states
		warn('You may want to use $exp->adaptor->store_states($exp)');
		$self->store_states($exp);
	  }
	}
	
    return \@args;
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
	$exp->mage_xml_id($sth->{'mysql_insertid'});

	$sql = "UPDATE experiment set mage_xml_id=".$exp->mage_xml_id()." where experiment_id =".$exp->dbID();
	$sth = $self->prepare($sql);
	$sth->execute();
  }

  return $exp;
}


=head2 list_dbIDs

  Args       : None
  Example    : my @exp_ids = @{$exp_a->list_dbIDs()};
  Description: Gets an array of internal IDs for all Experiment objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('experiment');
}


sub fetch_experiment_filter_counts{
  my $self = shift;

   my $sql = 'SELECT count(*), eg.name, eg.description, eg.is_project, ft.class, ct.name, ct.description '.
	'FROM experimental_group eg, experiment e, feature_set fs, feature_type ft, cell_type ct, '.
	  'status s, status_name sn '.
		'WHERE fs.experiment_id=e.experiment_id AND e.experimental_group_id=eg.experimental_group_id '.
		  'AND fs.feature_type_id=ft.feature_type_id AND fs.cell_type_id=ct.cell_type_id '.
			'AND fs.feature_set_id=s.table_id AND fs.type="annotated" AND s.table_name="feature_set" '.
			  'AND s.status_name_id=sn.status_name_id and sn.name="DISPLAYABLE" '.
				'GROUP BY eg.name, eg.is_project, ft.class, ct.name';


#  SELECT count(*), eg.name, eg.description, eg.is_project, ft.class, ct.name, ct.description FROM experimental_group eg, experiment e, feature_set fs, feature_type ft, cell_type ct, status s, status_name sn WHERE fs.experiment_id=e.experiment_id AND e.experimental_group_id=eg.experimental_group_id AND fs.feature_type_id=ft.feature_type_id AND fs.cell_type_id=ct.cell_type_id AND fs.feature_set_id=s.table_id AND fs.type="annotated" AND s.table_name="feature_set" AND s.status_name_id=sn.status_name_id and sn.name="DISPLAYABLE" and fs.name like "%NHEK%" GROUP BY eg.name, eg.is_project, ft.class, ct.name;
#This looks good, maybe the extra count are coming from below?


  my @rows = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
  my $ftype_info = $self->db->get_FeatureTypeAdaptor->get_regulatory_evidence_info;

  my %filter_info = ( 
					 #Project=> {},
					 #'Cell/Tissue' => {},
					 All =>
					 { All =>{ count       => 0,
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
	my $ft_class_label = $ftype_info->{$ft_class}{label};

	if(! exists $filter_info{'Evidence type'}{$ft_class_label}){
	  $filter_info{'Evidence type'}{$ft_class_label} = 
		{ count       => 0,
		  description => $ftype_info->{$ft_class}{long_name},
		};
	}
	$filter_info{'Evidence type'}{$ft_class_label}{count} += $count;
  }

  return \%filter_info;

  #Do we need to add an 'in_build' filter /data field?

}


1;

