#
# Ensembl module for Bio::EnsEMBL::Funcgen::OligoArray
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::OligoArray - A module to represent an oligonucleotide microarray.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::OligoArray;

my $array = Bio::EnsEMBL::Funcgen::OligoArray->new(
	-NAME        => 'Array-1',
        -FORMAT      => 'Tiled',
        -SIZE        => '1',
        -SPECIES     => 'Mus_musculus',
	-VENDOR      => 'Nimblegen',
        -DESCRIPTION => $desc,
);

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $array_adaptor = $db_adaptor->get_ArrayAdaptor();
my $array = $array_adaptor->fetch_by_name($array_name)

=head1 DESCRIPTION

An OligoArray object represents an oligonucleotide microarray. The data
(currently the name, format, size, species, vendor and description) are stored
in the array table.



=head1 AUTHOR

This module was created by Nathan Johnson, but is almost entirely based on the
AffyArray modules written by Ian Sealy and Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::OligoArray;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);# %VALID_TYPE);
@ISA = qw(Bio::EnsEMBL::Storable);


# Possible types for OligoArray objects
#%VALID_TYPE = (
#	'AFFY'  => 1,
#	'OLIGO' => 1,
#);


=head2 new

  Arg [-NAME]:
        string - the name of this array
  Arg [-TYPE]:
        string - the type of this array (AFFY or OLIGO)
  Example    : my $array = Bio::EnsEMBL::Funcgen::OligoArray->new(
								  -NAME        => 'Array-1',
								  -FORMAT      => 'Tiled',
								  -SIZE        => '1',
								  -SPECIES     => 'Mus_musculus',
								  -VENDOR      => 'Nimblegen',
								  -DESCRIPTION => $desc,
								 );
  Description: Creates a new Bio::EnsEMBL::Funcgen::OligoArray object.
  Returntype : Bio::EnsEMBL::Funcgen::OligoArray
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);

	#can we lc these?
	my ($name, $format, $size, $species, $vendor, $desc)
		= rearrange( ['NAME', 'FORMAT', 'SIZE', 'SPECIES', 'VENDOR', 'DESCRIPTION'], @_ );
	
	$self->name($name)         if defined $name;
	$self->format($format)     if defined $format;
	$self->size($size)         if defined $size;
	$self->species($species)   if defined $species;
	$self->vendor($vendor)     if defined $vendor;
	$self->description($desc)  if defined $desc;

	return $self;
}

=head2 get_all_Probes

  Args       : None
  Example    : my $probes = $array->get_all_Probes();
  Description: Returns all probes on an array. Needs a database connection.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoProbe objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_Probes {
	my $self = shift;

	if ( $self->dbID() && $self->adaptor() ) {
		my $opa = $self->adaptor()->db()->get_OligoProbeAdaptor();
		my $probes = $opa->fetch_all_by_Array($self);
		return $probes;
	} else {
		warning('Need database connection to retrieve Probes');
		return [];
	}
}

=head2 name

  Arg [1]    : (optional) string - the name of this array
  Example    : my $name = $array->name();
  Description: Getter, setter and lazy loader of name attribute for OligoArray
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub name {
	my $self = shift;
	$self->{'name'} = shift if @_;
	if ( !exists $self->{'name'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'name'};
}

=head2 format

  Arg [1]    : (optional) string - the format of the array
  Example    : my $format = $array->format();
  Description: Getter, setter and lazy loader of format attribute for
               OligoArray objects e.g. Tiled, Targetted etc...
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub format {
	my $self = shift;

	$self->{'format'} = shift if @_;
	if ( !exists $self->{'format'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'format'};
}

#is this going to be, # of array_chips, probe_sets or probes?
#currently number of array_chips...maybe change to probe_sets?
#Is array.size going to be used for validating array imports?
#Does it matter if we have an incomplete import from an experiment?
#No but we then need some way of checking further array imports to see what is missing

=head2 size

  Arg [1]    : (optional) int - the number of ? in the array
  Example    : my $size = $array->size();
  Description: Getter, setter and lazy loader of size attribute for
               OligoArray objects. The size is the number of ? in this array. 
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub size {
	my $self = shift;
	$self->{'size'} = shift if @_;
	if ( !exists $self->{'size'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'size'};
}

=head2 species

  Arg [1]    : (optional) string - the species of the array (e.g. Mus_musculus)
  Example    : my $species = $array->species();
  Description: Getter, setter and lazy loader of species attribute for OligoArray
               objects.
  Returntype : string
  Exceptions : Throws if argument cannot be mapped to a valid registry species alias
  Caller     : General
  Status     : Medium Risk

=cut

sub species {
  my $self = shift;
  my $species = shift;
  
  if ($species) {
    warn("Registry species alias lookup not yet implemented");
    $self->{'species'} = $species;
  }

  if ( !exists $self->{'species'} && $self->dbID() && $self->adaptor() ) {
    $self->adaptor->fetch_attributes($self);
  }
  return $self->{'species'};
}

=head2 vendor

  Arg [1]    : (optional) string - the name of the array vendor
  Example    : my $vendor = $array->vendor();
  Description: Getter, setter and lazy loader of vendor attribute for
               OligoArray objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub vendor {
	my $self = shift;
	$self->{'vendor'} = shift if @_;
	if ( !exists $self->{'vendor'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'vendor'};
}

=head2 description

  Arg [1]    : (optional) string - the description of the array
  Example    : my $size = $array->description();
  Description: Getter, setter and lazy loader of description attribute for
               OligoArray objects. 
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub description {
	my $self = shift;
	$self->{'description'} = shift if @_;
	if ( !exists $self->{'description'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'description'};
}



1;

