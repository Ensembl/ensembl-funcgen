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


my @true_tables = (['array', 'a']);
my @tables = @true_tables;

=head2 fetch_by_array_chip_dbID

  Arg [1]    : int - dbID of array_chip
  Example    : my $array = $oaa->fetch_by_array_chip_dbID($ac_dbid);
  Description: Retrieves a named Array object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

#Changed to use simple query extension
#Removed 1 query
#3.7 % or 1.04 times faster

sub fetch_by_array_chip_dbID {
  my ($self, $ac_dbid) = @_;
  
  throw('Must provide an ArrayChip dbID') if ! $ac_dbid;
  
  #Extend query tables
  push @tables, (['array_chip', 'ac']);
  
  #Extend query and group
  my $array = $self->generic_fetch('ac.array_chip_id='.$ac_dbid.' and ac.array_id=a.array_id GROUP by a.array_id')->[0];
  
  #Reset tables
  @tables = @true_tables;
  
  return $array;
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
    #my ($result) = @{$self->generic_fetch("a.name = '$name' and a.vendor='".uc($vendor)."'")};	
    return $self->generic_fetch("a.name = '$name' and a.vendor='".uc($vendor)."'")->[0];
}

=head2 fetch_by_name_class

  Arg [1]    : string - name of an array
  Arg [2]    : string - class of array e.g. AFFY_UTR
  Example    : my $array = $oaa->fetch_by_name_class('HuGene_1_0_st_v1', 'AFFY_ST');
  Description: Retrieves Array object from the database based on name and class.
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : Throws is name and class not passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name_class {
  my ($self, $name, $class) = @_;
  throw("Must provide and array and class e.g.'HuGene_1_0_st_v1', 'AFFY_ST'") if (! ($name && $class));
  #my ($result) = @{$self->generic_fetch("a.name = '$name' and a.class='".uc($class)."'")};	
  return $self->generic_fetch("a.name = '$name' and a.class='".uc($class)."'")->[0];
}


=head2 fetch_all_by_class

  Arg [1]    : string - class
  Example    : my $array = $oaa->fetch_all_by_class(''AFFY_ST');
  Description: Retrieves Array object from the database based class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : Throws if nor class passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_class {
    my ($self, $class) = @_;
    
    throw("Must provide and array class e.g.'AFFY_ST'") if (! defined $class);
    return $self->generic_fetch("a.class='".uc($class)."'");	
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

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my @arrays = @{$aa->fetch_all_by_Experiment($exp)};
  Description: Fetch all arrays associated with a given Experiment
               This is a convenience method to hide the 2 adaptors required 
               for this call.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : none
  Caller     : General
  Status     : at risk

=cut

#Changed to use simple query extension
#Removed 2 queries
#96.4% or 26.7 times faster!!!

sub fetch_all_by_Experiment{
  my ($self, $exp) = @_;

 $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exp);
  
  #Extend query tables
  push @tables, (['array_chip', 'ac'], ['experimental_chip', 'ec']);

  #Extend query and group
  my $arrays = $self->generic_fetch($exp->dbID.'=ec.experiment_id and ec.array_chip_id=ac.array_chip_id and ac.array_id=a.array_id GROUP by a.array_id');

  #Reset tables
  @tables = @true_tables;
	
  return $arrays;
}


=head2 fetch_all_by_ProbeSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::ProbeSet
  Example    : my @arrays = @{$aa->fetch_all_by_ProbeSet($probeset)};
  Description: Fetch all arrays containing a given ProbeSet
               This is a convenience method to hide the 2 adaptors required 
               for this call.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : none
  Caller     : General
  Status     : at risk

=cut


#Changed to use simple query extension
#Removed 1 query and hash loop
#This is only 1.04 times faster or ~ 4%

sub fetch_all_by_ProbeSet{
  my ($self, $pset) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ProbeSet', $pset);

  #Extend query tables
  push @tables, (['array_chip', 'ac'], ['probe', 'p']);

  #Extend query and group
  my $arrays =  $self->generic_fetch('p.probe_set_id='.$pset->dbID.' and p.array_chip_id=ac.array_chip_id and ac.array_id=a.array_id GROUP by a.array_id');

  #Reset tables
  @tables = @true_tables;
  
  return $arrays;
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
	
	return @tables;
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

