#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor - A database adaptor for fetching and
storing Funcgen Experiment objects.

=head1 SYNOPSIS

my $exp_a = $db->get_ExperimentAdaptor();
my $exp = $exp_a->fetch_by_name($name);


=head1 DESCRIPTION

The ExperimentAdaptor is a database adaptor for storing and retrieving
Funcgen Experiment objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning );
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_all_by_group

  Arg [1]    : int - dbID of array_chip
  Example    : my $array = $oaa->fetch_by_array_chip_dbID($ac_dbid);
  Description: Retrieves a named OligoArray object from the database.
  Returntype : listref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_group {
    my $self = shift;

    throw("Not yet implemented");

    my $ac_dbid = shift;
    my $sth = $self->prepare("
		SELECT a.array_id
		FROM array a, array_chip ac
		WHERE a.array_id = ac.array_id
		AND ac.array_chip_id = $ac_dbid
	");


    $sth->execute();
    my ($array_id) = $sth->fetchrow();
    
    return $self->fetch_by_dbID($array_id);
}



=head2 fetch_by_name

  Arg [1]    : string - name of an Experiment
  Example    : my $exp = $exp_a->fetch_by_name('Exp-1');
  Description: Retrieves a named Experiment object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : Throws if no name defined or if more than one returned
  Caller     : General
  Status     : At Risk -replace with fetch_all_by_name and fetch_by_name_group

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;
  
  throw("Need to specify and experiment name argument") if (! defined $name);

  my $result = $self->generic_fetch("e.name = '$name'");
  
  if (scalar @$result > 1) {
    throw("Experiment $name is not unique in the database, but only one result has been returned");
    #should have unique key of group_id and experiment_name
  } 
  return $result->[0];
}

=head2 get_all_experiment_names

  Arg [1]    : (optional) boolean - flag to denote whether experiment is flagged for web display
  Example    : my @names = $exp_a->get_all_experiment_names();
  Description: Retrieves names of all experiments.
  Returntype : ARRAYREF
  Exceptions : none
  Caller     : General
  Status     : At Risk -rename fetch?

=cut

sub get_all_experiment_names{
  my ($self, $displayable) = @_;

  my $sql = "SELECT e.name FROM experiment e";
  my @names = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  #can we do return [map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)}];
  return \@names;
}

#fetch_by_name_group


#=head2 fetch_all_by_design_type
#
#  Arg [1]    : List of strings - type(s) (e.g. AFFY or OLIGO)
#  Example    : my @arrays = @{$oaa->fetch_all_by_type('OLIGO')};
#  Description: Fetch all arrays of a particular type.
#  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoArray objects
#  Exceptions : Throws if no type is provided
#  Caller     : General
#  Status     : Medium Risk

#=cut

#sub fetch_all_by_design_type {
#	my ($self, @types) = @_;
	
#	throw('Need type as parameter') if !@types;
	
#	my $constraint;
#	if (scalar @types == 1) {
#		$constraint = qq( oa.type = '$types[0]' );
#	} else {
#		$constraint = join q(','), @types;
#		$constraint = qq( oa.type IN ('$constraint') );
#	}

#	return $self->generic_fetch($constraint);
#}

=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty Experiment objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Experiment getters
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
  Status     : At risk

=cut

sub _tables {
	my $self = shift;
	
	#should we add group, target, design_type, experimental_variable?

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
	
	return qw( e.experiment_id e.name e.egroup_id e.date e.primary_design_type e.description);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $exp_id, $name, $group_id, $p_design_type, $date, $description);
	
	$sth->bind_columns(\$exp_id, \$name, \$group_id, \$date, \$p_design_type, \$description);
	
	while ( $sth->fetch() ) {
	  my $exp = Bio::EnsEMBL::Funcgen::Experiment->new(
							   -DBID                => $exp_id,
							   -ADAPTOR             => $self,
							   -NAME                => $name,
							   -GROUP_ID            => $group_id,
							   -DATE                => $date,
							   -PRIMARY_DESIGN_TYPE => $p_design_type,
							   -DESCRIPTION         => $description,
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
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @args = @_;

	my ($s_exp);
   	
	my $sth = $self->prepare("INSERT INTO experiment
                                 (name, egroup_id, date, primary_design_type, description)
                                 VALUES (?, ?, ?, ?, ?)");
	

    foreach my $exp (@args) {
      if ( ! $exp->isa('Bio::EnsEMBL::Funcgen::Experiment') ) {
	warning('Can only store Experiment objects, skipping $exp');
	next;
      }
      
      my ($g_dbid) = $self->db->fetch_group_details($exp->group());
      throw("Group specified does, not exist.  Use Importer(group, location, contact)") if(! $g_dbid);
      
      # Has array already been stored?
      next if ( $exp->dbID() && $exp->adaptor() == $self );
      
      $s_exp = $self->fetch_by_name($exp->name());#validate on group too!
      
      if(! $s_exp){
	
	warn "Storing exp with name ".$exp->name()."\n";

	$sth->bind_param(1, $exp->name(),                SQL_VARCHAR);
	$sth->bind_param(2, $g_dbid,                    SQL_INTEGER);
	$sth->bind_param(3, $exp->date(),                SQL_VARCHAR);#date?
	$sth->bind_param(4, $exp->primary_design_type(), SQL_VARCHAR);
	$sth->bind_param(5, $exp->description(),         SQL_VARCHAR);
	
	$sth->execute();
	$exp->dbID($sth->{'mysql_insertid'});
	$exp->adaptor($self);
	
	
	#do we need to set egroup, target, design_type, experimentall_variable here?
      }
      else{
	warn("Experiment already exists in DB, using previously stored Experiment");
	$exp = $s_exp;
      }
    }
    
    return \@args;

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


#Need to add lazy load methods
#experimental_variables
#group


1;

