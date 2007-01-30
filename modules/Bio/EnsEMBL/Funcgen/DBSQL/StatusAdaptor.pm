#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::StatusAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::StatusAdaptor - A database adaptor for fetching
and setting the status of different entities in the database.

=head1 SYNOPSIS

my $status_a = $db->get_StatusAdaptor();

$status_a->set_status("DISPLAYABLE", $experiment);
my @displayable_states = @{$status_s->fetch_all_states_like("DISPLAYABLE", $experiment);

if($status_a->has_status($analysis_logic_name, $experiment){ ... }

if($status_a->has_status("DISPLAYABLE:".$analysis_logic_name, $experiment){ ... display methods ... }

my @exp_states = @{$status_a->fetch_all_states($experiment)};

=head1 DESCRIPTION

The StatusAdaptor is a database adaptor for fetching and setting the status 
of different entities in the database.  This provides a way to track whether an
experiment, experimental_chip or an array has been processed as: REGISTERED, IMPORTED etc...

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::StatusAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use vars qw(@ISA);
use strict;
use warnings;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 displayable_to_constraint

  Arg [1]    : Bio::EnsEMBL::Funcgen::DBSQL::"Adaptor"
  Arg [2]    : string (opt) - Constraint
  Arg [3]    : boolean - displayable
  Example    : $sql = $self->db->get_StatusAdaptor->displayable_to_constraint($self, $constraint, $displayable);
  Description: Appends the appropriate displayable constraint dependant on which adaptor is passed. 
  Returntype : string - constraint
  Exceptions : Throws if argument is not a Bio::EnsEMBL::Funcgen::DBSQL::"Adaptor"
  Caller     : Bio::EnsEMBL::Funcgen::DBSQL::"Adaptor" fetch methods
  Status     : At risk

=cut

sub displayable_to_constraint{
	my ($self, $adaptor, $constraint, $displayable) = @_;

	return $constraint if(! $displayable);

	if(! $adaptor || (ref($adaptor) !~ /Bio::EnsEMBL::Funcgen::DBSQL/)){
		throw("Need to specifiy a Bio::EnsEMBL::Funcgen::DBSQL::\"Adaptor\"");
	}

	my $table = $adaptor->_tables->[0];
	my $syn   = $adaptor->_tables->[1];
	
	my $sql = " AND ${syn}.${table}_id = status.table_id AND status.table_name='${table}' AND status.state='DISPLAYABLE'";

	return $sql;
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

  my $table = ${$obj->adaptor->_tables()}[0];

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

  my $sql = "SELECT state FROM status WHERE table_name=\"$table\" AND table_id=\"".$obj->dbID()."\"";

  my @states = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  return \@states;
}
             



=head2 has_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : Bio::EnsEMBL::"OBJECT"
  Example    : if($status_a->has_status('IMPORTED', $array){ ... skip import ... };
  Description: Tests wether a given object has a given state
  Returntype : BOOLEAN
  Exceptions : None
  Caller     : general
  Status     : At risk - somewhat defunct now due to Funcgen::Storable implementation

=cut



sub has_status{
  my ($self, $state, $obj) = @_;

  warn("To be deprecated, use Bio::EnsEMBL::Funcgen::Storable->has_status()");

  throw("cannot check state of an unstored object") if (! $obj->dbID());

  my $table = $self->_test_funcgen_table($obj);
  my $sql = "SELECT state FROM status WHERE table_name=\"$table\" AND table_id=\"".$obj->dbID()."\" AND state=\"$state\"";

  #could just return the call directly?
  my @row = $self->db->dbc->db_handle->selectrow_array($sql);

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

  if($self->has_status($state, $obj)){
    warning("$obj with dbID ".$obj->dbID()." already has state $state set\n");
  }else{
    my $table = $self->_test_funcgen_table($obj);
    my $sql = "INSERT INTO status(table_id, table_name, state) VALUES(\"".$obj->dbID()."\", \"$table\", \"$state\")";
    $self->db->dbc->do($sql);
  }

  return;
}

#quick method for ResultSetAdaptor->fetch_Resultfeatures_by_Slice

sub displayable_filter{
  my ($self, $table_name, @table_ids) = @_;

  throw("Must provide a table_name and table_ids to filter non-displayable ids") if(! ($table_name && @table_ids));
  
  my $sql = "SELECT table_id from status where table_name='$table_name' and table_id in (".join(", ", @table_ids).") and status.state='DISPLAYABLE'";
  
  
  my @displayable_ids = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  return \@displayable_ids;
	
}

1;

