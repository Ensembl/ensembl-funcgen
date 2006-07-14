#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor - A database adaptor for fetching and
storing Funcgen OligoArray objects.

=head1 SYNOPSIS

my $oaa = $db->get_OligoArrayAdaptor();

my $array = $oaa->fetch_by_name('Array-1');
my @arrays = @{$oaa->fetch_all()};

=head1 DESCRIPTION

The OligoArrayAdaptor is a database adaptor for storing and retrieving
Funcgen OligoArray objects.

=head1 AUTHOR

This module was created by Nathan Johnson, but is almost entirely based on the
OligoArrayAdaptor modules written by Ian Sealy and Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning );
use Bio::EnsEMBL::Funcgen::OligoArray;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_array_chip_dbID

  Arg [1]    : int - dbID of array_chip
  Example    : my $array = $oaa->fetch_by_array_chip_dbID($ac_dbid);
  Description: Retrieves a named OligoArray object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::OligoArray
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_array_chip_dbID {
    my $self = shift;
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

  Arg [1]    : string - name of an array
  Example    : my $array = $oaa->fetch_by_name('Array-1');
  Description: Retrieves a named OligoArray object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::OligoArray
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;
    
    my $result = $self->generic_fetch("a.name = '$name'");
	
    if (scalar @$result > 1) {
		warning("Array $name is not unique in the database, but only one result has been returned");
    } 


	#should have fetch by name vendor, to provide uniqueness?
	#should check for this on import!


    return $result->[0];
}

#=head2 fetch_all_by_type
#
#  Arg [1]    : List of strings - type(s) (e.g. AFFY or OLIGO)
#  Example    : my @arrays = @{$oaa->fetch_all_by_type('OLIGO')};
#  Description: Fetch all arrays of a particular type.
#  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoArray objects
#  Exceptions : Throws if no type is provided
#  Caller     : General
#  Status     : Medium Risk

#=cut

#sub fetch_all_by_type {
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

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoArray - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty OligoArray objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::OligoArray getters
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
  Status     : Medium Risk

=cut

sub _tables {
	my $self = shift;
	
	return ['array', 'a'];
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
	
	return qw( a.array_id a.name a.format a.size a.vendor a.description);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoArray objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $array_id, $name, $format, $size, $vendor, $description);
	
	$sth->bind_columns(\$array_id, \$name, \$format, \$size, \$vendor, \$description);
	
	while ( $sth->fetch() ) {
		my $array = Bio::EnsEMBL::Funcgen::OligoArray->new(
														   -dbID        => $array_id,
														   -adaptor     => $self,
														   -name        => $name,
														   -format      => $format,
														   -size        => $size,
														   -vendor      => $vendor,
														   -description => $description,
														  );

		push @result, $array;

		#if ($parent_id) {
		#	my $parent_array = Bio::EnsEMBL::Funcgen::OligoArray->new(
		#		-dbID    => $parent_id,
		#		-adaptor => $self,
		#	);
		#	$array->superset($parent_array);
		#}
	}
	return \@result;
}


=head2 _fetch_array_chips_by_array_dbID

  Arg [1]    : int - dbID or array
  Example    : None
  Description: Retrieves
  Returntype : Listref of array_chip hashes
  Exceptions : Throws if no array dbID specified
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _fetch_array_chips_by_array_dbID {
	my ($self, $array_dbID) = @_;

	throw("Must specifiy an array_id to retrieve array_chips") if (! $array_dbID);
	my (@result, $array_chip_id, $design_id, $name, %ac_tmp);

	my $sth = $self->prepare("select array_chip_id, design_id, name from array_chip where array_id = $array_dbID");	
	$sth->execute();

	$sth->bind_columns(\$array_chip_id, \$design_id, \$name);
	
	while ( $sth->fetch() ) {
		%ac_tmp = (
				   array_chip_id => $array_chip_id,
				   design_id => $design_id,
				   name => $name,
				  );

		push @result, {%ac_tmp};
	}
	return \@result;
}





=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::OligoArray objects
  Example    : $oaa->store($array1, $array2, $array3);
  Description: Stores given OligoArray objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub store {
    my $self = shift;
    my @args = @_;

	my $ac_sth = $self->prepare("INSERT INTO array_chip
                                 (design_id, array_id, name)
                                 VALUES (?, ?, ?)");


    foreach my $array (@args) {
		if ( !$array->isa('Bio::EnsEMBL::Funcgen::OligoArray') ) {
			warning('Can only store OligoArray objects, skipping $array');
			next;
		}

		if( ! exists $array->{'array_chips'}){
			warning("Need to implement array_chip inserts, should we create OligoArrayChip.pm? Skipping");
			next;
		}		

		# Has array already been stored?
		next if ( $array->dbID() && $array->adaptor() == $self );

		#my $superset = $array->superset();
		#if ( defined $superset && !$superset->dbID() ) {
		#	$self->store($superset);
		#}

		#can we prepare these statement once,and then bind in the loop, to prevent preparing multiple times?
		my $sth = $self->prepare("
			INSERT INTO array
			(name, format, size, vendor, description)
			VALUES (?, ?, ?, ?, ?)
		");
		$sth->bind_param(1, $array->name(),         SQL_VARCHAR);
		$sth->bind_param(2, $array->format(),       SQL_VARCHAR);
		$sth->bind_param(3, $array->size(),         SQL_INTEGER);
		$sth->bind_param(4, $array->vendor(),       SQL_VARCHAR);
		$sth->bind_param(5, $array->description(),  SQL_VARCHAR);

		$sth->execute();
		my $dbID = $sth->{'mysql_insertid'};
		$array->dbID($dbID);
		$array->adaptor($self);

		foreach my $array_chip(@{$array->array_chips()}){
			#validate keys here..or write object?!
			warn("Need to validate array_chip keys here or write obj");


			$ac_sth->bind_param(1, $array_chip->{'design_id'}, SQL_INTEGER);
			$ac_sth->bind_param(2, $dbID,                      SQL_INTEGER);
			$ac_sth->bind_param(3, $array_chip->{'name'},      SQL_VARCHAR);

			$ac_sth->execute();
			$$array_chip{'array_chip_id'} = $ac_sth->{'mysql_insertid'};

		}
	}
}

=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$oaa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoArray objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('array');
}


#New Funcgen methods
#fetch_all_by_group?
#fetch_by_channel_id?
#fetch_by_chip_id?
#fetch_by_probe_id
#fetch_by_probe_set_id


1;

