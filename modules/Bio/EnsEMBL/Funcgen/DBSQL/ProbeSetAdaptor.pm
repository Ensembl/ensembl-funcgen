#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::ProbeSetAdaptor - A database adaptor for fetching and
storing ProbeSet objects.

=head1 SYNOPSIS

my $opa = $db->get_ProbeSetAdaptor();

my $probeset = $opa->fetch_by_array_probeset_name('Array-1', 'ProbeSet-1');

=head1 DESCRIPTION

The ProbeSetAdaptor is a database adaptor for storing and retrieving
ProbeSet objects.

=head1 AUTHOR

This module was created by Nathan Johnson, but is almost entirely based on the
ProbeAdaptor module written by Ian Sealy and Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Probe;#?
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


#may need to pass array object, as there is a possibilty of it being non-unique between vendors?

=head2 fetch_by_array_probeset_name

  Arg [1]    : string - name of array
  Arg [2]    : string - name of probeset
  Example    : my $probeset = $opsa->fetch_by_array_probeset_name('Array-1', 'Probeset-1');
  Description: Returns a probeset given the array name and probeset name
               This will uniquely define a probeset. Only one
			   probeset is ever returned.
  Returntype : Bio::EnsEMBL::ProbeSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_array_probeset_name {
	my $self          = shift;
	my $array_name    = shift;
	my $probeset_name = shift;
	
	my $sth = $self->prepare("
		SELECT probe_set_id
		FROM probe_set ps, array a, array_chip ac
		WHERE a.array_id = ac.array_id
		AND a.name = ?
		AND ps.name = ?
	");
	
	$sth->bind_param(1, $array_name,    SQL_VARCHAR);
	$sth->bind_param(2, $probeset_name, SQL_VARCHAR);

	$sth->execute();
	
	my ($probeset_id) = $sth->fetchrow();
	
	if ($probeset_id) {
		return $self->fetch_by_dbID($probeset_id);
	} else {
		return undef;
	}
}



=head2 fetch_all_by_Array

  Arg [1]    : Bio::EnsEMBL::Array
  Example    : my @probesets = @{$opsa->fetch_all_by_Array($array)};
  Description: Fetch all probes on a particular array.
  Returntype : Listref of Bio::EnsEMBL::ProbeSet objects.
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Array {
	my $self  = shift;
	my $array = shift;

	throw("Not yet implemented");

	my ($probeset_id, @probesets);
	
	if ( !ref($array) || !$array->isa('Bio::EnsEMBL::Array') ) {
		warning('fetch_all_by_Array requires a Bio::EnsEMBL::Array object');
		return [];
	}
	
	my $array_id = $array->dbID();
	if (!defined $array_id) {
		warning('fetch_all_by_Array requires a stored Bio::EnsEMBL::Array object');
		return [];
	}


	#Nath
	#retrieve all array_chip_ids and do a generic fetch using a joined or statement?
	#or
	#build and array of probesets using the fetch_by_dbID method

	
	my $sth = $self->prepare("
		SELECT probe_set_id
		FROM probe_set ps, array a, array_chip ac
		WHERE a.array_id = ac.array_id
        AND ac.array_chip_id = ps.array_chip_id
		AND a.name = $array_id
	");
	

	$sth->execute();
	
	while($probeset_id = $sth->fetchrow()){
		push @probesets, $self->fetch_by_dbID($probeset_id);
	}

	return \@probesets;
}

=head2 fetch_by_ProbeFeature

  Arg [1]    : Bio::EnsEMBL::ProbeFeature
  Example    : my $probeset = $opsa->fetch_by_ProbeFeature($feature);
  Description: Returns the probeset that created a particular feature.
  Returntype : Bio::EnsEMBL::ProbeSet
  Exceptions : Throws if argument is not a Bio::EnsEMBL::ProbeFeature object
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_ProbeFeature {
	my $self    = shift;
	my $feature = shift;
	
	if (
		!ref($feature)
		|| !$feature->isa('Bio::EnsEMBL::Funcgen::ProbeFeature')
		|| !$feature->{'probe_id'}
	) {
		throw('fetch_by_ProbeFeature requires a stored Bio::EnsEMBL::Funcgen::ProbeFeature object');
	}
	
	my $sth = $self->prepare("
		SELECT probe_set_id
		FROM probe_set ps, probe p, probe_feature pf
		WHERE pf.probe_id = p.probe_id
        AND ps.probe_set_id = p.probe_set_id
		AND pf.probe_feature_id = ?
	");

	$sth->bind_param(1, $feature->{'probe_id'},    SQL_VARCHAR);
	$sth->execute();

	my ($probeset_id) = $sth->fetchrow();

	return $self->fetch_by_dbID($probeset_id);
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

  	return [ 'probe_set', 'ps' ];
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

  #remove xref_id and use xref tables
  return qw( ps.probe_set_id ps.name ps.size ps.family);

}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates ProbeSet objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::ProbeSet objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@result, $current_dbid, $probeset_id, $name, $size, $family);
	my ($array, %array_cache);
	
	$sth->bind_columns( \$probeset_id,  \$name, \$size, \$family);
	

	#do not have array_chip adaptor
	#use array adaptor directly
	#how are we going ot handle the cache here?????

	my $probeset;
	while ( $sth->fetch() ) {
		#$array = $array_cache{$array_id} || $self->db->get_ArrayAdaptor()->fetch_by_dbID($array_id);

		#This is nesting array object in probeset!
		#$array = $array_cache{$arraychip_id} || $self->db->get_ArrayAdaptor()->fetch_by_array_chip_dbID($arraychip_id);

		#Is this required? or should we lazy load this?
		#Should we also do the same for probe i.e. nest or lazy load probeset
		#Setting here prevents, multiple queries, but if we store the array cache in the adaptor we can overcome this
		#danger of eating memory here, but it's onld the same as would be used for generating all the probesets
		#what about clearing the cache?
		#also as multiple array_chips map to same array, cache would be redundant
		#need to store only once and reference.
		#have array_cache and arraychip_map
		#arraychip_map would give array_id which would be key in array cache
		#This is kinda reinventing the wheel, but reducing queries and redundancy of global cache
		#cache would never be populated if method not called
		#there for reducing calls and memory, increasing speed of generation/initation
		#if method were called
		#would slightly slow down processing, and would slightly increase memory as cache(small as non-redundant)
		#and map hashes would persist

		#Do we even need this????

		warn("Can we lazy load the arrays from a global cache, which is itself lazy loaded and non-redundant?\n");

		
		#this current id stuff is due to lack of probeset table in core
		#if (!$current_dbid || $current_dbid != $probeset_id) {
		  
		  # New probeset
		  $probeset = Bio::EnsEMBL::Funcgen::ProbeSet->new
			(
			 -dbID         => $probeset_id,
			 -name         => $name,
			 -size         => $size,
			 #			 -array        => $array,
			 -family       => $family,
			 -adaptor     => $self,
			);
		push @result, $probeset;

			#$current_dbid = $probeset_id;
		#} else {
		#	# Extend existing probe
		#	$probe->add_Array_probename($array, $name);
		#}
	}
	return \@result;
}

=head2 store

  Arg [1]    : List of Bio::EnsEMBL::Funcgen::ProbeSet objects
  Example    : $opa->store($probeset1, $probeset2, $probeset3);
  Description: Stores given ProbeSet objects in the database. Should only be
               called once per probe because no checks are made for duplicates.??? It certainly looks like there is :/
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if arguments are not Probe objects
  Caller     : General
  Status     : At Risk

=cut

sub store {
	my ($self, @probesets) = @_;

	my ($sth, $array);

	if (scalar @probesets == 0) {
		throw('Must call store with a list of Probe objects');
	}

	my $db = $self->db();

	PROBESET: foreach my $probeset (@probesets) {

		if ( !ref $probeset || !$probeset->isa('Bio::EnsEMBL::Funcgen::ProbeSet') ) {
			throw('ProbeSet must be an ProbeSet object');
		}

		if ( $probeset->is_stored($db) ) {
			warning('ProbeSet [' . $probeset->dbID() . '] is already stored in the database');
			next PROBESET;
		}
		
		# Get all the arrays this probe is on and check they're all in the database
		#my $arrays = $probeset->get_all_Arrays();
		#my @stored_arrays;
		#for $array (@$arrays) {
		#	if ( defined $array->dbID() ) {
		#		push @stored_arrays, $array;
		#	}
		#}
		#if ( !@stored_arrays ) {
		#	warning('ProbeSets need attached arrays to be stored in the database');
		#	next PROBESET;
		#}

		# Insert separate entry (with same probe_set_id) in probe_set
		# for each array the probe is on
		#my $dbID;
		#for $array (@stored_arrays) {
		#		

		#	if (defined $dbID) {
		#		$sth = $self->prepare("
		#			INSERT INTO probe_set
		#			(probe_set_id, name, size, array_chip_id, family)
		#			VALUES (?, ?, ?, ?, ?, ?)
		#		");
		#

		#		# Probesets we've seen already
		#		$sth->bind_param(1, $dbID,                                 SQL_INTEGER);
		#		$sth->bind_param(2, $probeset->name(),                     SQL_VARCHAR);
		#		$sth->bind_param(3, $probeset->size(),                     SQL_INTEGER);
		#		$sth->bind_param(4, $probeset->array_chip_id(),            SQL_INTEGER);
		#		$sth->bind_param(5, $probeset->family(),                   SQL_VARCHAR);
		#		$sth->execute();
		#	} else {
			
		# New probeset
		$sth = $self->prepare("
					INSERT INTO probe_set
					(name, size, family)
					VALUES (?, ?, ?)
				");
		$sth->bind_param(1, $probeset->name(),                     SQL_VARCHAR);
		$sth->bind_param(2, $probeset->size(),                     SQL_INTEGER);
		$sth->bind_param(3, $probeset->family(),                   SQL_VARCHAR);
		
		$sth->execute();
		my $dbID = $sth->{'mysql_insertid'};
		$probeset->dbID($dbID);
		$probeset->adaptor($self);
		#}
		#}
	}

	return \@probesets;
}

=head2 list_dbIDs

  Arg [1]    : none
  Example    : my @ps_ids = @{$opa->list_dbIDs()};
  Description: Gets an array of internal IDs for all ProbeSet objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my ($self) = @_;

	return $self->_list_dbIDs('probe_set');
}

1;

