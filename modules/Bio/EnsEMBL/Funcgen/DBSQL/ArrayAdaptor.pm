#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor - A database adaptor for fetching and
storing Funcgen Array objects.

=head1 SYNOPSIS

my $oaa = $db->get_ArrayAdaptor();

my $array = $oaa->fetch_by_name('Array-1');
my @arrays = @{$oaa->fetch_all()};

=head1 DESCRIPTION

The ArrayAdaptor is a database adaptor for storing and retrieving
Funcgen Array objects.

=head1 AUTHOR

This module was created by Nathan Johnson, but is almost entirely based on the
ArrayAdaptor modules written by Ian Sealy and Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

=head2 fetch_by_array_chip_dbID

  Arg [1]    : int - dbID of array_chip
  Example    : my $array = $oaa->fetch_by_array_chip_dbID($ac_dbid);
  Description: Retrieves a named Array object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : None
  Caller     : General
  Status     : At Risk

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
  Description: Retrieves a named Array object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;


    throw("This method is deprecated, use fetch_by_name_vendor");
    
    my $result = $self->generic_fetch("a.name = '$name'");
	
    if (scalar @$result > 1) {
		warning("Array $name is not unique in the database, but only one result has been returned");
    } 


    #should have fetch by name vendor, to provide uniqueness?
    #should check for this on import!
    return $result->[0];
}



=head2 fetch_by_name_vendor

  Arg [1]    : string - name of an array
  Arg [2]    : string - name of vendor e.g. NIMBLEGEN
  Example    : my $array = $oaa->fetch_by_name('Array-1');
  Description: Retrieves a named Array object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name_vendor {
    my ($self, $name, $vendor) = @_;
    
    throw("Must provide and array and vendor name") if (! ($name && $vendor));

    #unique key means this only returns one element
    my ($result) = @{$self->generic_fetch("a.name = '$name' and a.vendor='".uc($vendor)."'")};	
    return $result;
}



=head2 fetch_all_by_type

  Arg [1]    : List of strings - type(s) (e.g. OLIGO, PCR)
  Example    : my @arrays = @{$aa->fetch_all_by_type('OLIGO')};
  Description: Fetch all arrays of a particular type.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : Throws if no type is provided
  Caller     : General
  Status     : at risk

=cut

sub fetch_all_by_type {
  my ($self, @types) = @_;
	
  throw('Need type as parameter') if ! @types;
	
  my $constraint;
  if (scalar @types == 1) {
    $constraint = qq( a.type = '$types[0]' );
  } else {
    $constraint = join q(','), @types;
    $constraint = qq( a.type IN ('$constraint') );
  }

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Experiment

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiement
  Example    : my @arrays = @{$aa->fetch_all_by_Experiment($exp)};
  Description: Fetch all arrays associated with a given Experiment
               This is a convenience method to hide the 3 adaptor required 
               for this call.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : none
  Caller     : General
  Status     : at risk

=cut

sub fetch_all_by_Experiment{
  my ($self, $exp) = @_;

  my %array_ids;
	
  my $echips = $self->db->get_ExperimentalChipAdaptor->fetch_all_by_Experiment($exp);

  foreach my $achip(@{$self->db->get_ArrayChipAdaptor->fetch_all_by_ExperimentalChips($echips)}){
	$array_ids{$achip->array_id()} = 1;
  }

  return $self->generic_fetch('a.array_id IN ('.join(', ', keys %array_ids).')');
}


=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Funcgen::Array - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty Array objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Array getters
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
  Status     : At Risk

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
  Status     : At Risk

=cut

sub _columns {
	my $self = shift;
	
	return qw( a.array_id a.name a.format a.vendor a.description a.type a.class);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;
	
  my (@result, $array_id, $name, $format, $vendor, $description, $type, $class);
  
  $sth->bind_columns(\$array_id, \$name, \$format, \$vendor, \$description, \$type, \$class);
  
  while ( $sth->fetch() ) {

    my $array = Bio::EnsEMBL::Funcgen::Array->new
	  (
	   -dbID        => $array_id,
	   -adaptor     => $self,
	   -name        => $name,
	   -format      => $format,
	   -vendor      => $vendor,
	   -description => $description,
	   -type        => $type,
	   -class       => $class, 
	  );

    push @result, $array;
    
 
  }
  return \@result;
}




=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::Array objects
  Example    : $oaa->store($array1, $array2, $array3);
  Description: Stores given Array objects in the database. This
               method checks for arrays previously stored and updates 
               and new array_chips accordingly.
  Returntype : Listref of stored Array objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


#This works slightly differently as arary_chip are not stored as objects, 
#yet we need to retrieve a dbID for the array before we know about all the array_chips

sub store {
  my $self = shift;
  my @args = @_;
  
  my ($sarray);
  
  my $sth = $self->prepare("
			INSERT INTO array
			(name, format, vendor, description, type, class)
			VALUES (?, ?, ?, ?, ?, ?)");
  
    
  foreach my $array (@args) {
    if ( !$array->isa('Bio::EnsEMBL::Funcgen::Array') ) {
      warning('Can only store Array objects, skipping $array');
      next;
    }
    
    if (!( $array->dbID() && $array->adaptor() == $self )){
      #try and fetch array here and set to array if valid
      $sarray = $self->fetch_by_name_vendor($array->name(), $array->vendor());#this should be name_vendor?
      
      if( ! $sarray){
		#sanity check here
		throw("Array name must not be longer than 30 characters") if (length($array->name) > 40);
		$sth->bind_param(1, $array->name(),         SQL_VARCHAR);
		$sth->bind_param(2, $array->format(),       SQL_VARCHAR);
		$sth->bind_param(3, $array->vendor(),       SQL_VARCHAR);
		$sth->bind_param(4, $array->description(),  SQL_VARCHAR);
		$sth->bind_param(5, $array->type(),         SQL_VARCHAR);
		$sth->bind_param(6, $array->class(),        SQL_VARCHAR);
		
		
		$sth->execute();
		my $dbID = $sth->{'mysql_insertid'};
		$array->dbID($dbID);
		$array->adaptor($self);
      }
      else{
		#warn("Array already stored, using previously stored array\n");# validating array_chips\n");
		$array = $sarray;
      }
    }
  }
    
  return \@args;
}




=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$oaa->list_dbIDs()};
  Description: Gets an array of internal IDs for all Array objects in the
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

