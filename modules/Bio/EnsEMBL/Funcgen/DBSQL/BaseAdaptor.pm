
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

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use DBI qw(:sql_types);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor Exporter);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});


#do we want to keep the IMPORTED_CS status for feature_sets/array_chips?
#rename to MAPPED_CS_N?


=head2 store_states

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->store_states($result_set);
  Description: Stores states of Storable in status table.
  Returntype : None
  Exceptions : Throws if Storable is not stored
  Caller     : General
  Status     : At Risk

=cut


sub store_states{
  my ($self, $storable) = @_;

  throw('Must pass a Bio::EnsEMBL::Funcgen::Storable') if(! $storable->isa("Bio::EnsEMBL::Funcgen::Storable"));
  
  foreach my $state(@{$storable->get_all_states()}){

    $self->store_status($state, $storable) if (! $self->has_stored_status($state, $storable));
  }

  return;

}

=head2 fetch_all

  Arg[1]     : string - optional status name e.g. 'DISPLAYABLE'
  Example    : my @dsets = @{$dsa->fetch_all()};
  Description: Gets all available objects from the DB, which 
               might not be a good idea, shouldnt be called on 
               the BIG tables though
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all {
  my ($self, $status) = @_;

  my ($constraint);
  #Can we throw here if we're trying to get all from known large tables


  $constraint = $self->status_to_constraint($status) if $status;

  
  if(defined $status && ! defined $constraint){
	#warn "You have specifed a status($status) which is not present in the DB";
	return undef;
  }


  return $self->generic_fetch($constraint);

}



=head2 fetch_all_diplayable

  Example    : my @displayable_dset = @{$dsa->fetch_all_displayable()};
  Description: Gets all displayable DataSets
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk - can we just reimplement fetch_all with optional status arg

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
  Status     : At Risk - To be removed

=cut

sub fetch_all_by_status{ 
  my ($self, $status) = @_; 

  deprecate('Use fetch_all($status) instead');
  return $self->fetch_all($status);

}


#Can we not just re implement fetch_all here to have the optional status arg?



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
  my $status_id = $self->_get_status_name_id($status);

  
  #THIS DOES NOT ACCOMODATE THE EXPEIRMENTAL_SUBSET ISSUE!!

  #NO we need to handle this better
  #can't just return as we'd then simply ignore the contraint
  #can't throw??
  
  return if (! $status_id);
  
  my @tables = $self->_tables;
  my ($table_name, $syn) = @{$tables[0]};

  my @status_ids;

  my $sql = "SELECT table_id from status where table_name='$table_name' and status_name_id='$status_id'";
  @status_ids = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  #This is causing problems as we might get none, which will invalidate the sql
  #Hence we return nothing

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
  #ExpreimetnalSubSet fix, as doesn't have own adaptor
  $table = 'experimental_subset' if $obj->isa('Bio::EnsEMBL::Funcgen::ExperimentalSubset');


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
  my $status_id = $self->_get_status_name_id($state);

  if (! (ref($obj) && $obj->isa("Bio::EnsEMBL::Funcgen::Storable") && $obj->dbID())){
	throw("Must pass a stored Bio::EnsEMBL::Funcgen::Storable");
  }

  my $table = $self->_test_funcgen_table($obj);
 


  if($status_id){
    my $sql = "SELECT status_name_id FROM status WHERE table_name=\"$table\" AND table_id=\"".$obj->dbID()."\" AND status_name_id=\"$status_id\"";

    #could just return the call directly?
    @row = $self->db->dbc->db_handle->selectrow_array($sql);
  }

  return (@row) ? 1 : 0;
}




=head2 store_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : Bio::EnsEMBL::"OBJECT"
  Example    : $status_a->store_status('IMPORTED', $array_chip);
  Description: Sets a state for a given object
  Returntype : None
  Exceptions : Warns if state already set
               Throws is status name is not already stored.
  Caller     : general
  Status     : At risk - Move to Status

=cut


sub store_status{
  my ($self, $state, $obj) = @_;

  my $sql;

  if($self->has_stored_status($state, $obj)){
    warn("$obj with dbID ".$obj->dbID()." already has status $state set\n");
  }else{
    my $status_id = $self->_get_status_name_id($state);

    if(! $status_id){
      throw("$state is not a valid status_name for $obj:\t".$obj->dbID);
      #$sql = "INSERT into status_name(name) values('$state')";
      #$self->db->dbc->do($sql);
      #$status_id = $self->_get_status_name_id($state);
    }

    my $table = $self->_test_funcgen_table($obj);
  

    $sql = "INSERT into status(table_id, table_name, status_name_id) VALUES('".$obj->dbID()."', '$table', '$status_id')";
    $self->db->dbc->do($sql);

	#Should we not be setting it in the obj here too?
	#No becasue we should have already added to the object.
  }

  return;
}


=head2 revoke_status

  Arg [1]    : string - status name e.g. 'IMPORTED'
  Arg [2]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->revoke_status('DAS DISPLAYABLE', $result_set);
  Description: Revokes the given state of Storable in status table.
  Returntype : None
  Exceptions : Warns if storable does not have state
               Throws is status name is not valid
               Throws if not state passed
               Throws if Storable is not valid and stored
  Caller     : General
  Status     : At Risk

=cut


sub revoke_status{
  my ($self, $state, $storable) = @_;

  throw('Must provide a status name') if(! defined $state);
  throw('Must pass a Bio::EnsEMBL::Funcgen::Storable') if(! $storable->isa("Bio::EnsEMBL::Funcgen::Storable"));
 
  my $status_id = $self->_get_status_name_id($state);
  my $table_name = $self->_test_funcgen_table($storable);
  #hardcode for ExperimentalSubset as this uses the ExperimentalSetAdaptor
  $table_name = 'experimental_subset' if $storable->isa('Bio::Ensembl::Funcgen:ExperimentalSubset');
 

  if(! $self->has_store_status($state, $storable)){
	warn $storable.' '.$storable->dbID()." does not have status $state to revoke\n";
	return;
  }

  #do sanity checks on table to ensure that IMPORTED does not get revoke before data deleted?
  #how would we test this easily?

  my $sql = "delete from status where table_name='${table_name}'".
	" and status_name_id=${status_id} and table_id=".$storable->dbID();

  $self->db->dbc->db_handle->do($sql);

  #now splice from status array;
  #splice in loop should work as we will only see 1
  #Just hash this?

  for my $i(0..$#{$storable->{'states'}}){
	
	if($storable->{'states'}->[0] eq $state){
	  splice @{$storable->{'states'}}, $i, 1;
	  last;
	}
  }

  return;
}

=head2 revoke_states

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->revoke_status($result_set);
  Description: Revokes all states of Storable in status table.
  Returntype : Bio::EnsEMBL::Funcgen::Storable
  Exceptions : None
  Caller     : General + Helper rollback methods
  Status     : At Risk

=cut


sub revoke_states{
  my ($self, $storable) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Storable', $storable);
 
  my $table_name = $self->_test_funcgen_table($storable);

  #hardcode for ExperimentalSubset as this uses the ExperimentalSetAdaptor
  $table_name = 'experimental_subset' if $storable->isa('Bio::Ensembl::Funcgen:ExperimentalSubset');
 
  my $sql = "delete from status where table_name='${table_name}'".
	" and table_id=".$storable->dbID();

  $self->db->dbc->db_handle->do($sql);

  #Clear stored states
  undef $storable->{'states'};

  return $storable;
}



=head2 status_filter

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Arg [2]    : string - table name e.g. experimental_chip
  Arg [3]    : list   - table dbIDs
  Exmaple    : my @displayable_ec_ids = @{$ec_adaptor->status_filter('DISPLAYABLE', 
                                                                     'experimental_chip', 
                                                                     (map $_->dbID, @echips))};
  Description: Quick method for filtering dbIDs based on their table and and status
  Returntype : ARRAYREF
  Exceptions : Warns if state already set
               Throws is status name is not already stored.
  Caller     : general - ResultSetAdaptor->fetch_Resultfeatures_by_Slice
  Status     : At risk - Move to Status?

=cut


sub status_filter{
  my ($self, $status, $table_name, @table_ids) = @_;


  my @status_ids;

  my $status_id = $self->_get_status_name_id($status);


  return \@status_ids if(! $status_id);

  throw("Must provide a table_name and table_ids to filter non-displayable ids") if(! ($table_name && @table_ids));
  
  my $sql = "SELECT table_id from status where table_name='$table_name' and table_id in (".join(", ", @table_ids).") and status.status_name_id='$status_id'";
  
  
  @status_ids = map $_ = "@$_", @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  return \@status_ids;
	
}


=head2 _get_status_name_id

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Example    : my $status_id = $self->_get_status_name_id('IMPORTED');
  Description: Retrieves the dbID of a given status_name
  Returntype : INT
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::BaseAdaptor
  Status     : At risk - Move to Status?

=cut


sub _get_status_name_id{
  my ($self, $status) = @_;

  #$self->_validate_status($status);

  my $sql = "SELECT status_name_id from status_name where name='$status'";

  my $ref = $self->db->dbc->db_handle->selectrow_arrayref($sql);

  my ($status_id) = @$ref if $ref;


  #we should throw here?
  #To force manual addition of the status_name
  #need to make sure all status_names which are explicitly used by API
  #are stored in all DBs, else we could find ourselves with broken code
  #for sparsely populated DBs

  throw("Status name $status is not valid.  Maybe you need to add it using update_status_name.pl?") if ! $status_id;


  return $status_id;
}




=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               An external identifier of the feature to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Example    : my @features =
                  @{ $adaptor->fetch_all_by_external_name( 'NP_065811.1') };
  Description: Retrieves all features which are associated with
               an external identifier such as a GO term, Swissprot
               identifer, etc.  Usually there will only be a single
               feature returned in the list reference, but not
               always.  Features are returned in their native
               coordinate system, i.e. the coordinate system in which
               they are stored in the database.  If they are required
               in another coordinate system the Feature::transfer or
               Feature::transform method can be used to convert them.
               If no features with the external identifier are found,
               a reference to an empty list is returned.
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

#This might be mor eefficient if we wrote DBEntryAdaptor->fetch_all_by_

sub fetch_all_by_external_name {
  my ( $self, $external_name, $external_db_name ) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();
  my (@ids, $type, $type_name);
  ($type = ref($self)) =~ s/.*:://;
  $type =~ s/Adaptor$//;
  ($type_name = $type) =~ s/Feature$/_feature/;
  my $xref_method = 'list_'.lc($type_name).'_ids_by_extid';

  if(! $entryAdaptor->can($xref_method)){
	warn "Does not yet accomodate $type external names";
	return;
  }
  else{
	@ids = $entryAdaptor->$xref_method($external_name, $external_db_name);
  }

  return $self->fetch_all_by_dbID_list( \@ids );
}





1;

