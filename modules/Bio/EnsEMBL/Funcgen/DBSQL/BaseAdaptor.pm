
#
# Perl module for Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor - Simple wrapper class for Funcgen Storable Adaptors

=head1 SYNOPSIS

$adaptor->store_states($storable);

=head1 DESCRIPTION

This is a simple wrapper class to hold common methods to all Funcgen StorableAdaptors.
Includes status methods.

=head1 SEE ALSO

Bio::EnsEMBL::DBSQL::BaseAdaptor

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



=head2 compose_constraint_query

  Arg [1]    : Hash - Params hash containing 'constraints' key value pairs
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_FeatureType($type);
  Description: Retrieves FeatureSet objects from the database based on feature_type id.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if arg is not a valid FeatureType
  Caller     : General
  Status     : At Risk

=cut

#This approach originated in the need for a fully flexible method call
#for the FeatureSetAdaptor to support the Experiment view.
#Add support for final_clause?
#Will still be useful to define an array of valid constraints in the descendant adaptor
#This will enable us to restrict the generic constraints(in here) for a given adaptor
#and dynamically provide a list of valid constraint in the error output

sub compose_constraint_query{
  my ($self, $params) = @_;

  if($params &&
	 (ref($params) ne 'HASH') ){
	throw('You must pass a valid params HASHREF to compose_constraint_query');
  }
  
  
  my @constraints;

  if( exists ${$params}{constraints}){
		
	my @filter_names = keys (%{$params->{constraints}});

	foreach my $constraint_key(keys (%{$params->{constraints}})){

    my $constrain_method = '_constrain_'.$constraint_key;

    if(! $self->can($constrain_method)){
      throw($constraint_key." is not a valid constraint");

      #Need to add test on and list valid constraints

      # please specify values for one of:\t".
      #           join(', ', keys(%constraint_config)));
    }

    my ($constraint, $constraint_conf) = $self->$constrain_method($params->{constraints}{$constraint_key});
    push @constraints, $constraint;
    

    #Currently only handle tables here but could also 
    #set other dynamic config e.g. final clause

    if (exists ${$constraint_conf}{tables}) {
      push @{$self->TABLES}, @{$constraint_conf->{tables}};
    }
	} # END OF CONSTRAINTS
  }	# END OF $PARAMS				

  return join(' AND ', @constraints) || '';
}



sub reset_true_tables{
  my $self = shift;

  #deref to avoid modifying TRUE_TABLES
  @{$self->TABLES} = @{$self->TRUE_TABLES};
  return;
}

=head2 store_states

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : $rset_adaptor->store_states($result_set);
  Description: Stores states of Storable in status table.
  Returntype : None
  Exceptions : Throws if Storable is not stored
  Caller     : General
  Status     : At Risk

=cut

#do we want to keep the IMPORTED_CS status for feature_sets/array_chips?
#rename to MAPPED_CS_N?

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
	my @tables = $self->_tables;
	my ($table_name, undef) = @{$tables[0]};

	warn "You are trying to fetch $table_name entries which have a $status status, which is not present in the DB";
	return [];
  }

  return $self->generic_fetch($constraint);

}



=head2 fetch_all_displayable

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
  
  #This is now supported in compose_constraint_query
  #Which avoid some of the problems below

  my $constraint;

  #This will throw if status not valid, but still may be absent
  my $status_id = $self->_get_status_name_id($status);

    
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

  #Can we change this to is_stored_and_valid for a Storable?
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Storable', $obj);
  #Does this test for ad

  my @tables = $obj->adaptor->_tables;

  my ($table) = @{$tables[0]};
  #InputSubSet fix, as doesn't have own adaptor
  $table = 'input_subset' if $obj->isa('Bio::EnsEMBL::Funcgen::InputSubset');

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
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::BaseAdaptor
  Status     : At risk

=cut



sub has_stored_status{
  my ($self, $state, $obj) = @_;

  my (@row);

  #Only used for set_status, merge with set_status?
  my $status_id = $self->_get_status_name_id($state);
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
  Exceptions : None
  Caller     : general
  Status     : At risk - Move to Status.pm?

=cut


sub store_status{
  my ($self, $state, $obj) = @_;

  my $sql;
  my $table = $self->_test_funcgen_table($obj);

  if(! $self->has_stored_status($state, $obj)){
    my $status_id = $self->_get_status_name_id($state);

    if(! $status_id){
      throw("$state is not a valid status_name for $obj:\t".$obj->dbID);
    }

	$sql = "INSERT into status(table_id, table_name, status_name_id) VALUES('".$obj->dbID()."', '$table', '$status_id')";
	$self->db->dbc->do($sql);

	#Setting it in the obj if it is not already present.
	$obj->add_status($state) if(! $obj->has_status($state, $obj));
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
  Caller     : General
  Status     : At Risk

=cut


sub revoke_status{
  my ($self, $state, $storable) = @_;

  throw('Must provide a status name') if(! defined $state);
  my $table_name = $self->_test_funcgen_table($storable);
  my $status_id = $self->_get_status_name_id($state);
  
  #hardcode for InputSubset
  $table_name = 'input_subset' if $storable->isa('Bio::Ensembl::Funcgen:InputSubset');
 

  if(! $self->has_stored_status($state, $storable)){
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


=head2 set_imported_states_by_Set

  Arg [1]    : Bio::EnsEMBL::Funcgen::Set e.g. a FeatureSet or ResultSet
  Example    : $self->set_imported_states_by_Set($set);
  Description: Sets default states for imported Feature|ResultSets
  Returntype : None
  Exceptions : None
  Caller     : Import parsers and RunnableDBs
  Status     : At risk - move to BaseImporter

=cut

#This needs to be used by RunnableDBs too!
#All state stuff is handled by BaseAdaptor?
#Can we put this in the SetAdaptor?

sub set_imported_states_by_Set{
  my ($self, $set) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Set', $set);
  #This should really be restricted to FeatureSet and ResultSet

  #Store default states for FeatureSets
  #DAS_DISPLAYABLE IMPORTED_'CSVERSION'
  #These need to insert ignore as may already be present?
  #Insert ignore may not catch an invalid status
  #So add states and store states as this checks
  $set->adaptor->store_status('DAS_DISPLAYABLE', $set);
  

  #To get assembly version here we need to
  # 1 get the current default chromosome version
  # or/and 
  # 2 Use the assembly param to guess it from the coord_sys table
  # #This may pose problems for DB names which use numbers in their genebuild version
  # Would need to set this as a flag and/or specify the genebuild version too
  # Currently dnadb is set to last dnadb with 'assembly' as default version
  # We should match as test, just to make sure
  
  #Get default chromosome version for this dnadb
  my $cs_version = $self->db->dnadb->get_CoordSystemAdaptor->fetch_by_name('chromosome')->version;

  #Sanity check that it matches the assembly param?
  #Woould need to do this if ever we loaded on a non-default cs version

  $set->adaptor->store_status("IMPORTED_${cs_version}", $set);
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

  throw("Status name $status is not valid.  Maybe you need to add it to the status_name table?") if ! $status_id;


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
               an external identifier such as an Ensembl Gene or Transcript
               stable ID etc.  Usually there will only be a single
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
  Exceptions : Warns if method not available for given object adaptor
  Caller     : general
  Status     : at risk

=cut

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
	return [];
  }
 
  #Would be better if _list_ method returned and arrayref
  return $self->fetch_all_by_dbID_list([$entryAdaptor->$xref_method($external_name, $external_db_name)]);
}


=head2 fetch_all_by_external_names

  Arg [1]    : ARRAYREF of strings. External identifiers of the features to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Example    : my @features =
                  @{ $adaptor->fetch_all_by_external_names(['ENST00003548913', ...])};
  Description: Retrieves all features which are associated with
               the external identifiers such as a Ensembl gene or transcript 
               stable IDs, etc.  Features are returned in their native
               coordinate system, i.e. the coordinate system in which
               they are stored in the database.  If they are required
               in another coordinate system the Feature::transfer or
               Feature::transform method can be used to convert them.
               If no features with the external identifier are found,
               a reference to an empty list is returned.
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : Warns if xref method not available for given object adaptor
  Caller     : general
  Status     : at risk

=cut

sub fetch_all_by_external_names{
  my ( $self, $external_names, $external_db_name ) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();
  my ($type, $type_name);
  ($type = ref($self)) =~ s/.*:://;
  $type =~ s/Adaptor$//;
  ($type_name = $type) =~ s/Feature$/_feature/;
  my $xref_method = 'list_'.lc($type_name).'_ids_by_extids';
  

  if(! $entryAdaptor->can($xref_method)){
	warn "Does not yet accomodate $type external names";
	return [];
  }
  
 
  #Would be better if _list_ method returned and arrayref
  my @ids = $entryAdaptor->$xref_method($external_names, $external_db_name);

  return $self->fetch_all_by_dbID_list(\@ids);
}


=head2 fetch_all_by_linked_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript
  Example    : my @psets =
                  @{ $probe_set_adaptor->fetch_all_by_linked_Transcript($tx_obj);
  Description: Retrieves all features which are associated with
               the given Ensembl Transcript.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : Throws if arguments not valid
  Caller     : general
  Status     : at risk

=cut

sub fetch_all_by_linked_Transcript{
  my ($self, $tx) = @_;

  if(! $tx ||
	 ! (ref($tx) && $tx->isa('Bio::EnsEMBL::Transcript') && $tx->dbID)){
	throw('You must provide a valid stored Bio::EnsEMBL:Transcript object');
  }
  
  return $self->fetch_all_by_external_name($tx->stable_id, $self->db->species.'_core_Transcript')
}

=head2 fetch_all_by_linked_transcript_Gene

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : my @psets =
                  @{ $probe_set_adaptor->fetch_all_by_linked_transcript_Gene($gene_obj);
  Description: Retrieves all features which are indirectly associated with
               the given Ensembl Gene, through it's Transcripts.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Storable objects
               Maybe any Feature, FeatureType, Probe or ProbeSet
  Exceptions : Throws if arguments not valid
  Caller     : general
  Status     : at risk

=cut


sub fetch_all_by_linked_transcript_Gene{
   my ( $self, $gene ) = @_;

   if(! $gene ||
	  ! (ref($gene) && $gene->isa('Bio::EnsEMBL::Gene') && $gene->dbID)){
	 throw('You must provide a valid stored Bio::EnsEMBL:Gene object');
   }
   #No need to quote param here as this is a known int from the DB.
   my $tx_sids = $gene->adaptor->db->dbc->db_handle->selectcol_arrayref('select tsid.stable_id from transcript_stable_id tsid, transcript t where t.gene_id='.$gene->dbID.' and t.transcript_id=tsid.transcript_id');
   
   return $self->fetch_all_by_external_names($tx_sids, $self->db->species.'_core_Transcript');
}



=head2 store_associated_feature_types

  Arg [1]     : Bio::EnsEMBL::Funcgen::Sotrable
  Example     : $ext_feat_adaptor->store_associated_feature_type($ext_feat);
  Description : Store FeatureTypes assoicated with a given Storable
  Returntype  : None
  Exceptions  : Throws if FeatureTypes are not valid or stored
  Caller      : Adaptors
  Status      : At risk

=cut 

sub store_associated_feature_types { 
  my ($self, $storable) = @_;

  #Direct access to avoid lazy loading with an unstored SetFeature
  my $assoc_ftypes = $storable->{'associated_feature_types'};

  #Could be undef or empty
  return if ! defined $assoc_ftypes || scalar(@$assoc_ftypes) == 0;

  my $table_name = $storable->adaptor->_main_table->[0];
  my $dbid = $storable->dbID;

  my $sql = 'INSERT into associated_feature_type(table_id, table_name, feature_type_id) values (?,?,?)';

  foreach my $ftype(@$assoc_ftypes){

	#We have already tested the class but not whether it is stored
	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);

	my $sth = $self->prepare($sql);
	$sth->bind_param(1, $dbid,        SQL_INTEGER);
	$sth->bind_param(2, $table_name,  SQL_VARCHAR);
	$sth->bind_param(3, $ftype->dbID, SQL_INTEGER);
	$sth->execute();
  }	


  return;
} 


=head2 fetch_all_by_associated_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : my $assoc_ftypes = $ft_adaptor->fetch_all_by_associated_SetFeature($ext_feature);
  Description: Fetches all objects which have associated FeatureType. 
               Note this is not the main FeatureType for this object.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Storable objects
  Exceptions : Throws if FeatureType not valid or stored
  Caller     : General
  Status     : At risk

=cut


sub fetch_all_by_associated_FeatureType{
  my ($self, $ftype) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);
  my ($table_name, $table_syn) = @{$self->_main_table};

  push @{$self->TABLES}, ['associated_feature_type', 'aft'];
  my $constraint = "aft.feature_type_id=? AND aft.table_name='${table_name}' AND aft.table_id=${table_syn}.${table_name}_id";

  $self->bind_param_generic_fetch($ftype->dbID,  SQL_INTEGER);
  my $objs = $self->generic_fetch($constraint);
  $self->reset_true_tables;

  return $objs;
}


=head2 _main_table

  Example    : my $syn = $adaptor->_main_table->[1];
  Description: Convenience method to retrieve the main table or main table synonym for this adaptor
               Entirely dependent on ensembl convention of always having main table as first element
               of tables array.
  Returntype : Array ref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub _main_table{
  my $self = shift;

  #Need to do this to put it in list context to avoid just returning the last value
  my @tables = $self->_tables();
  return $tables[0];
}


=head2 _list_dbIDs

  Example    : my @table_ids = @{$adaptor->_list_dbIDs()};
  Description: Wrapper for parent method, dynamically passes correct table name to query.
               Gets an array of internal IDs for all objects in the main table of this class.
  Returntype : List of Ints
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut


sub _list_dbIDs{
  my $self = shift;
  return $self->SUPER::_list_dbIDs($self->_main_table->[0]);
}



#Need this to support generating valid adaptor names in FeatureSet
#this needs to be in the BaseAdaptor, so we can access it without 
#having a set or feature adaptor defined.


sub build_feature_class_name{
  my ($self, $fclass) = @_;
   
  if(! defined $fclass){
    throw('You must pass a feature class argument to build the feature class name');
  }
 
  my @words = split('_', $fclass);
  
  foreach my $word(@words){
    $word = ucfirst($word);
    $word = 'DNA' if $word eq 'Dna';
  }
  
  return join('', (@words, 'Feature') );
}



# Generic constraint methods, called from compose_constraint_query
# Depends on inheriting Adaptor having TABLES and TRUE_TABLES 'constants' set.

sub _constrain_status {
  my ($self, $state) = @_;
  
  my $constraint_conf = { tables => [['status', 's']]};  #,['status_name', 'sn']),
    
  #This can't use IN without duplicating the result
  #Would need to add a default_final_clause to group

  #my @sn_ids;
  #if( (ref($states) ne 'ARRAY') ||
  #scalar(@$states) == 0 ){
  #throw('Must pass an arrayref of status_names');
  #}
  #foreach my $sn(@$states){
  ##This will throw if status not valid, but still may be absent
  #	push @sn_ids, $self->_get_status_name_id($sn);
  #  }
      
	  
      
  my @tables = $self->_tables;
  my ($table_name, $syn) = @{$tables[0]};
	  
  my $constraint = " $syn.${table_name}_id=s.table_id AND ".
    "s.table_name='$table_name' AND s.status_name_id=".$self->_get_status_name_id($state);
  #"s.table_name='$table_name' AND s.status_name_id IN (".join(', ', @sn_ids.')';
  
  return ($constraint, $constraint_conf);
}

### DEPRECATED ###

sub list_dbIDs { #Deprecated in v69
	my $self = shift;	
  deprecate('Please use _list_dbIDs.');
	return $self->_list_dbIDs($self->_main_table->[0]);
}



1;

