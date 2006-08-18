#
# Ensembl module for Bio::EnsEMBL::DBSQL::OligoProbeAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::OligoProbeAdaptor - A database adaptor for fetching and
storing OligoProbe objects.

=head1 SYNOPSIS

my $opa = $db->get_OligoProbeAdaptor();

my $probe = $opa->fetch_by_array_probe_probeset('Array-1', 'Probe-1', undef);

=head1 DESCRIPTION

The OligoProbeAdaptor is a database adaptor for storing and retrieving
OligoProbe objects.

=head1 AUTHOR

This module was created by Ian Sealy, but is almost entirely based on the
OligoProbeAdaptor module written by Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::OligoProbeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::OligoProbe;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_array_probe_probeset_name

  Arg [1]    : string - name of array
  Arg [2]    : string - name of probe
  Arg [3]    : (optional) string - name of probeset
  Example    : my $probe = $opa->fetch_by_array_probeset_probe('Array-1', 'Probe-1', 'ProbeSet-1');
  Description: Returns a probe given a combination of array name, probeset and
               probe name. This will uniquely define an Affy probe. Only one
			   probe is ever returned.
  Returntype : Bio::EnsEMBL::Funcgen::OligoProbe
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_array_probe_probeset_name {
	my $self          = shift;
	my $array_name    = shift;
	my $probe_name    = shift;
	my $probeset_name = shift;
	
	my $probeset_clause = "";
	my $ops_table_alias = "";

	if(defined $probeset_name){
		$probeset_clause = "AND (op.oligo_probe_set_id = ops.oligo_probe_set_id AND ops.name = $probeset_name)";
		$ops_table_alias = ", oligo_probe_set ops";
	}

	# Need to deal with non-Affy probes where probeset is NULL
	# (Also allow probeset to be empty string, just in case)
	#if (!$probeset_name) {
	#	$probeset_name = '';
	#	$probeset_clause = "(op.oligo_probe_set_id IS NULL OR $probeset_clause)";
	#}
	
	#need to do a look up of all ac_ids and do a joined OR statement
	my $array_ref = $self->dbc->db_handle->selectall_arrayref("select ac.array_chip_id from array_chip ac, array a where a.array_id = ac.array_id and a.name = \"$array_name\"");
	map $_ = "@{$_}", @$array_ref;#only works for one element arrays, as we're really turning it into a space separated string.
	my $ac_clause = "(op.array_chip_id = ".join(" OR op.array_chip_id = ", @$array_ref).")";

	my $sth = $self->prepare("
		SELECT op.oligo_probe_id
		FROM oligo_probe op $ops_table_alias
		WHERE $ac_clause
		$probeset_clause
		AND op.name = ?
	");
	
	#$sth->bind_param(1, $array_name,    SQL_VARCHAR);
	#$sth->bind_param(1, $probeset_name, SQL_VARCHAR);
	$sth->bind_param(1, $probe_name,    SQL_VARCHAR);
	$sth->execute();
	
	my ($probe_id) = $sth->fetchrow();

	if ($probe_id) {
		return $self->fetch_by_dbID($probe_id);
	} else {
		return undef;
	}
}

=head2 fetch_all_by_probeset

  Arg [1]    : string - probeset name
  Example    : my @probes = @{$opa->fetch_all_by_probeset('Probeset-1')};
  Description: Fetch all probes in a particular probeset.
  Returntype : Listref of Bio::EnsEMBL::OligoProbe objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_probeset {
	my $self     = shift;
	my $probeset = shift;

	#use Oligo_ProbeSet adaptor?
	#my $probe_set_id = $self->db->db_handle->selectrow_array("select oligo_probe_set_id from oligo_probe_set);
	my $probe_set_id = $self->db->get_OligoProbeSetAdaptor->fetch_by_name($probeset)->probe_set_id();
	return $self->generic_fetch("op.oligo_probe_set_id = '$probe_set_id'");
}

=head2 fetch_all_by_Array

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoArray
  Example    : my @probes = @{$opa->fetch_all_by_Array($array)};
  Description: Fetch all probes on a particular array.
  Returntype : Listref of Bio::EnsEMBL::OligoProbe objects.
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_Array {
	my $self  = shift;
	my $array = shift;
	
	if ( !ref($array) || !$array->isa('Bio::EnsEMBL::Funcgen::OligoArray') ) {
		warning('fetch_all_by_Array requires a Bio::EnsEMBL::Funcgen::OligoArray object');
		return [];
	}
	
	my $array_id = $array->dbID();
	if (!defined $array_id) {
		warning('fetch_all_by_Array requires a stored Bio::EnsEMBL::Funcgen::OligoArray object');
		return [];
	}

	#get all array_chip_ids, for array and do a multiple OR statement with generic fetch
	
	return $self->generic_fetch("op.array_chip_id = ".join(" OR op.array_chip_id = ", @{$array->get_array_chip_ids()}));
}

=head2 fetch_by_OligoFeature

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoFeature
  Example    : my $probe = $opa->fetch_by_OligoFeature($feature);
  Description: Returns the probe that created a particular feature.
  Returntype : Bio::EnsEMBL::OligoProbe
  Exceptions : Throws if argument is not a Bio::EnsEMBL::Funcgen::OligoFeature object
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_OligoFeature {
	my $self    = shift;
	my $feature = shift;
	
	if (
		!ref($feature)
		|| !$feature->isa('Bio::EnsEMBL::Funcgen::OligoFeature')
		|| !$feature->{'_probe_id'}
	) {
		throw('fetch_by_OligoFeature requires a stored Bio::EnsEMBL::Funcgen::OligoFeature object');
	}
	
	return $self->fetch_by_dbID($feature->{'_probe_id'});
}

=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _tables {
	my $self = shift;

  	return [ 'oligo_probe', 'op' ];
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _columns {
  my $self = shift;

  return qw( op.oligo_probe_id op.oligo_probe_set_id op.name op.length op.array_chip_id op.class);

}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoProbe objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoProbe objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $current_dbid, $arraychip_id, $probe_id, $array_id, $probe_set_id, $name, $class, $probelength);
	my ($array, %array_cache, %probe_set_cache);
	
	$sth->bind_columns(\$probe_id, \$probe_set_id, \$name, \$probelength, \$arraychip_id, \$class);
	
	my $probe;
	while ( $sth->fetch() ) {

		warn("Need to sort array cacheing, have redundant cache!!");
		#This is nesting array and probeset objects in probe!
		

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


		$array = $array_cache{$arraychip_id} || $self->db->get_OligoArrayAdaptor()->fetch_by_array_chip_dbID($arraychip_id);

		
		#I don't think we need this?  Certainly not for storing

		#$probe_set = $probe_set_cache{$probe_set_id} || $self->db->get_OligoArrayAdaptor()->fetch_by_array_chip_dbID($arraychip_id);
		#probe_set cache would be substantially bigger!!
		#potentially as many as the probes

		#Just build cache and nest for now,may want to just return ID and lazy load

		my ($probeset);

		if($probe_set_id){
			$probeset = $probe_set_cache{$probe_set_id} || $self->db->get_OligoProbeSetAdaptor()->fetch_by_dbID($probe_set_id);
		}

		if (!$current_dbid || $current_dbid != $probe_id) {
			# New probe

			#UC??? or does rearrange handle this?

			$probe = Bio::EnsEMBL::Funcgen::OligoProbe->new
			  (
			   -dbID          => $probe_id,
			   -name          => $name,
			   -array_chip_id => $arraychip_id,
			   -array         => $array,
			   -probe_set     => $probeset,
			   -length        => $probelength,
			   -class         => $class,
			   -adaptor       => $self,
			);
			push @result, $probe;
			$current_dbid = $probe_id;
		} else {
			# Extend existing probe
			$probe->add_array_chip_probename($arraychip_id, $name);
		}
	}
	return \@result;
}

=head2 store

  Arg [1]    : List of Bio::EnsEMBL::Funcgen::OligoProbe objects
  Example    : $opa->store($probe1, $probe2, $probe3);
  Description: Stores given OligoProbe objects in the database. Should only be
               called once per probe because no checks are made for duplicates
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if arguments aren't OligoProbe objects
  Caller     : General
  Status     : Medium Risk

=cut

sub store {
	my ($self, @probes) = @_;

	my ($ac_id, $sth);

	if (scalar @probes == 0) {
		throw('Must call store with a list of OligoProbe objects');
	}

	my $db = $self->db();

	PROBE: foreach my $probe (@probes) {

		if ( !ref $probe || !$probe->isa('Bio::EnsEMBL::Funcgen::OligoProbe') ) {
			throw('Probe must be an OligoProbe object');
		}

		if ( $probe->is_stored($db) ) {
			warning('OligoProbe [' . $probe->dbID() . '] is already stored in the database');
			next PROBE;
		}
		
		# Get all the arrays this probe is on and check they're all in the database
		my %array_hashes;

		foreach $ac_id (keys %{$probe->{'arrays'}}) {

			if (defined ${$probe->{'arrays'}}{$ac_id}->dbID()) {
				$array_hashes{$ac_id} = $probe->{'arrays'}{$ac_id};
			}
		}

		if ( ! %array_hashes ) {
			warning('Probes need attached arrays to be stored in the database');
			next PROBE;
		}

		# Insert separate entry (with same oligo_probe_id) in oligo_probe
		# for each array/array_chip the probe is on
		my $dbID;

		foreach $ac_id (keys %array_hashes) {			
			my $ps_id = (defined $probe->probeset()) ? $probe->probeset()->dbID() : undef;

			if (defined $dbID) {
				# Probe we've seen already
				$sth = $self->prepare("
					INSERT INTO oligo_probe
					( oligo_probe_id, oligo_probe_set_id, name, length, array_chip_id, class )
					VALUES (?, ?, ?, ?, ?, ?)
				");
				
				$sth->bind_param(1, $dbID,                                                     SQL_INTEGER);
				$sth->bind_param(2, $ps_id,                                                    SQL_INTEGER);
				$sth->bind_param(3, $probe->get_probename($array_hashes{$ac_id}->name()),    SQL_VARCHAR);
				$sth->bind_param(4, $probe->length(),                                          SQL_INTEGER);
				$sth->bind_param(5, $ac_id,                                                    SQL_INTEGER);
				$sth->bind_param(6, $probe->class(),                                           SQL_VARCHAR);
				$sth->execute();
			} else {
				# New probe
				$sth = $self->prepare("
					INSERT INTO oligo_probe
					( oligo_probe_set_id, name, length, array_chip_id, class )
					VALUES (?, ?, ?, ?, ?)
				");

				$sth->bind_param(1, $ps_id,                                                    SQL_INTEGER);
				$sth->bind_param(2, $probe->get_probename($array_hashes{$ac_id}->name()),    SQL_VARCHAR);
				$sth->bind_param(3, $probe->length(),                                          SQL_INTEGER);
				$sth->bind_param(4, $ac_id,                                                    SQL_INTEGER);
				$sth->bind_param(5, $probe->class(),                                           SQL_VARCHAR);

			
				$sth->execute();
				$dbID = $sth->{'mysql_insertid'};
				$probe->dbID($dbID);
				$probe->adaptor($self);
			}
		}
	}
}

=head2 list_dbIDs

  Arg [1]    : none
  Example    : my @feature_ids = @{$opa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoProbe objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my ($self) = @_;

	return $self->_list_dbIDs('oligo_probe');
}

1;

