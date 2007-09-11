
#
# BioPerl module for Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor
##
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor - Simple wrapper class for Funcgen StorableAdaptors

=head1 SYNOPSIS

$adaptor->store_states($storable);

=head1 DESCRIPTION

This is a simple wrapper class to hold common methods to all Funcgen StorableAdaptors.
Includes status methods.

Post questions to the EnsEMBL developer mailing list: <ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;
require Exporter;
use vars qw(@ISA @EXPORT);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use DBI qw(:sql_types);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor Exporter);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});




=head2 store_states

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->store_states($result_set);
  Description: Stores states of Storable in status table.
  Returntype : None
  Exceptions : Throws if Storable is not stored
  Caller     : General
  Status     : At Risk

=cut


#We need to control what can be stored, so we need to check cell_line_ids?
#Or is this level of control to be implicit?  Would there be a use for a multi-celled DataSet
#Yes!  So we let the DB take anything, and make the obj_from_sth method handle all


sub store_states{
  my ($self, $storable) = @_;

  throw('Must call store with a list of OligoFeature objects') if(! $storable->isa("Bio::EnsEMBL::Funcgen::Storable"));
  
  foreach my $state(@{$storable->get_all_states()}){

    $self->set_status($state, $storable) if (! $self->has_stored_status($state, $storable));
  }

  return;

}


=head2 fetch_all_diplayable

  Example    : my @displayable_dset = @{$dsa->fetch_all_displayable()};
  Description: Gets all displayable DataSets
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_displayable{
  my $self = shift;
  return $self->fetch_all_by_status('DISPLAYABLE');
}

=head2 fetch_all_by_status

  Arg [1]    : string - status e.g. 'DISPLAYABLE'
  Example    : my @displayable_dset = @{$dsa->fetch_all_by_status('DISPLAYABLE')};
  Description: Gets all DataSets with given status
  Returntype : ARRAYREF
  Exceptions : Throws is no status defined
               Warns if  
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_status{ 
  my ($self, $status) = @_; 
  my $constraint = $self->status_to_constraint('DISPLAYABLE');

  return (defined $constraint) ? $self->generic_fetch($constraint) : undef;
}


=head2 status_to_constraint

  Arg [1]    : string - status e.g. 'DISPLAYABLE'
  Arg [2]    : string - Constraint
  Example    : $sql = $self->status_to_constraint($self, $constraint, $status);
  Description: Appends the appropriate status constraint dependant on the BaseAdaptor sub class.
  Returntype : string - constraint
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::DBSQL::"BaseAdaptors"
  Status     : At risk

=cut

sub status_to_constraint{
  my ($self, $status) = @_;
  
  my $constraint;

  #This will throw if status not valid, but still may be absent
  my $status_id = $self->get_status_id($status);

  


  #NO we need to handle this better
  #can't just return as we'd then simply ignore the contraint
  #can't throw??
  

  return if (! $status_id);
  
  my @tables = $self->_tables;
  my ($table_name, $syn) = @{$tables[0]};

  my @status_ids;

  my $sql = "SELECT table_id from status where table_name='$table_name' and status_name_id='$status_id'";
    
  @status_ids = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
  $constraint = " $syn.${table_name}_id IN (".join(',', @status_ids).")" if @status_ids;
  

  return $constraint;

}



=head2 _test_funcgen_table

  Arg [1]    : Bio::EnsEMBL::"OBJECT"
  Example    : $status_a->_is_funcgen_object($experimental_chip)};
  Description: Tests if the object is a valid funcgen object with an identifiable table_name
  Returntype : string - table_name
  Exceptions : Throws if argument if a Bio::EnsEMBL::Funcgen::"OBJECT" not supplied
               Throws if not table name identified
  Caller     : general
  Status     : At risk

=cut



sub _test_funcgen_table{
  my ($self, $obj) = @_;

  if(! $obj || 
     (ref($obj) !~ /Bio::EnsEMBL::Funcgen/) ||
     (ref($obj) =~ /::DBSQL::/)){
    throw("Need to specifiy a Bio::EnsEMBL::Funcgen::\"OBJECT\"");
  }

  throw("Cannot test state of unstored object: $obj") if (! $obj->is_stored($self->db()));

  my @tables = $self->_tables;

  my ($table) = @{$tables[0]};

  #my $table = ${$obj->adaptor->_tables()}[0];

  return $table || $self->throw("Cannot identify table name from $obj adaptor");
}



=head2 fetch_all_states

  Arg [1]    : Bio::EnsEMBL::"OBJECT"
  Arg [2...] : listref of states
  Example    : my @ec_states = @{$status_a->fetch_all_states($experimental_chip)};
  Description: Retrieves all states associated with the given "OBJECT"
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut

sub fetch_all_states{
  my ($self, $obj) = @_;

  my $table = $self->_test_funcgen_table($obj);

  my $sql = "SELECT name FROM status_name sn, status s WHERE s.table_name='$table' AND s.table_id='".$obj->dbID()."' and s.status_name_id=sn.status_name_id";


  my @states = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  return \@states;
}
             



=head2 has_stored_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : Bio::EnsEMBL::Storable
  Example    : if($status_a->has_stored_status('IMPORTED', $array){ ... skip import ... };
  Description: Tests wether a given object has a given state
  Returntype : BOOLEAN
  Exceptions : Throws if Storable not passed or stored
  Caller     : Bio::EnsEMBL::Funcgen::BaseAdaptor
  Status     : At risk

=cut



sub has_stored_status{
  my ($self, $state, $obj) = @_;

  my (@row);

  #Only used for set_status, merge with set_status?
  my $status_id = $self->get_status_id($state);

  throw("Must pass a stored Bio::EnsEMBL::Funcgen::Storable") if (! ($obj->isa("Bio::EnsEMBL::Funcgen::Storable") && $obj->dbID()));

  my $table = $self->_test_funcgen_table($obj);


  if($status_id){
    my $sql = "SELECT status_name_id FROM status WHERE table_name=\"$table\" AND table_id=\"".$obj->dbID()."\" AND status_name_id=\"$status_id\"";

    #could just return the call directly?
    @row = $self->db->dbc->db_handle->selectrow_array($sql);
  }

  return (@row) ? 1 : 0;
}

=head2 set_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : Bio::EnsEMBL::"OBJECT"
  Example    : $status_a->set_status('IMPORTED', $array_chip);
  Description: Sets a state for a given object
  Returntype : None
  Exceptions : Warns if state already set
  Caller     : general
  Status     : At risk - Move to Status

=cut


sub set_status{
  my ($self, $state, $obj) = @_;

  my $sql;

  if($self->has_stored_status($state, $obj)){
    warn("$obj with dbID ".$obj->dbID()." already has status $state set\n");
  }else{
    my $status_id = $self->get_status_id($state);

    if(! $status_id){
      warn("Creating NEW status_name entry for $state.  Is this a valid state?");
      $sql = "INSERT into status_name(name) values('$state')";
      $self->db->dbc->do($sql);
      $status_id = $self->get_status_id($state);
    }

    my $table = $self->_test_funcgen_table($obj);
    
    $sql = "INSERT into status(table_id, table_name, status_name_id) VALUES('".$obj->dbID()."', '$table', '$status_id')";
    $self->db->dbc->do($sql);
  }

  return;
}

#quick method for ResultSetAdaptor->fetch_Resultfeatures_by_Slice

sub status_filter{
  my ($self, $status, $table_name, @table_ids) = @_;


  my @status_ids;

  my $status_id = $self->get_status_id($status);


  return \@status_ids if(! $status_id);

  throw("Must provide a table_name and table_ids to filter non-displayable ids") if(! ($table_name && @table_ids));
  
  my $sql = "SELECT table_id from status where table_name='$table_name' and table_id in (".join(", ", @table_ids).") and status.status_name_id='$status_id'";
  
  
  @status_ids = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  return \@status_ids;
	
}

sub get_status_id{
  my ($self, $status) = @_;

  $self->_validate_status($status);

  my $sql = "SELECT status_name_id from status_name where name='$status'";

  my $ref = $self->db->dbc->db_handle->selectrow_arrayref($sql);

  my ($status_id) = @$ref if $ref;

  return $status_id;
}


sub _validate_status{
  my ($self, $status) = @_;

  throw("Must pass a status to validate") if ! $status;

  my $valid = 0;

  #We could do some look up on the table here, but this may compound problems if someone has hacked the table

  my @state_regexs = ('IMPORTED', 'IMPORTED_CS_', 'DISPLAYABLE', 'RESOLVED');


  foreach my $regex(@state_regexs){
    $valid = 1 if ($status =~ /$regex/);
  }

  throw("Not a valid status: $status") if(! $valid);
  
  return;
}

1;

