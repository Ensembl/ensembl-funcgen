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
	-INCLUDED_IN => $another_array,
	-SETSIZE     => 1,
	-TYPE        => 'OLIGO',
);

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

package Bio::EnsEMBL::OligoArray;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA %VALID_TYPE);
@ISA = qw(Bio::EnsEMBL::Storable);


# Possible types for OligoArray objects
%VALID_TYPE = (
	'AFFY'  => 1,
	'OLIGO' => 1,
);


=head2 new

  Arg [-NAME]:
        string - the name of this array
  Arg [-INCLUDED_IN]:
        (optional) Bio::EnsEMBL::OligoArray - a possible superset array
  Arg [-SETSIZE]: 
        int - the number of probes in a probe set (usually 1 unless Affy)
  Arg [-TYPE]: 
        string - the type of this array (AFFY or OLIGO)
  Example    : my $array = Bio::EnsEMBL::OligoArray->new(
                   -NAME        => 'Array-1',
				   -INCLUDED_IN => $another_array,
				   -SETSIZE     => 1,
				   -TYPE        => 'OLIGO',
               );
  Description: Creates a new Bio::EnsEMBL::OligoArray object.
  Returntype : Bio::EnsEMBL::OligoArray
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);

	my ($name, $superset, $setsize, $type)
		= rearrange( ['NAME', 'INCLUDED_IN', 'SETSIZE', 'TYPE'], @_ );
	
	$self->name($name)         if defined $name;
	$self->superset($superset) if defined $superset;
	$self->setsize($setsize)   if defined $setsize;
	$self->type($type)         if defined $type;

	return $self;
}

=head2 get_all_Probes

  Args       : None
  Example    : my $probes = $array->get_all_Probes();
  Description: Returns all probes on an array. Needs a database connection.
  Returntype : Listref of Bio::EnsEMBL::OligoProbe objects
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

=head2 superset

  Arg [1]    : (optional) Bio::EnsEMBL::OligoArray - a superset array
  Example    : my $parent_array = $array->superset();
  Description: Getter, setter and lazy loader of superset attribute for
               OligoArray objects. A superset is another OligoArray that
               contains all the probes of this OligoArray. This is bordering
               on superfluous.
  Returntype : Bio::EnsEMBL::OligoArray
  Exceptions : Throws if argument isn't a Bio::EnsEMBL::OligoArray object
  Caller     : General
  Status     : Medium Risk

=cut

sub superset {
	my $self = shift;
	my $superset = shift;
	if ($superset) {
		if (!ref $superset || !$superset->isa('Bio::EnsEMBL::OligoArray')) {
			throw('Superset must be a Bio::EnsEMBL::OligoArray');
		}
		$self->{'superset'} = $superset;
	}
	if ( !exists $self->{'superset'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'superset'};
}

=head2 setsize

  Arg [1]    : (optional) int - the number of probes in a probe set
  Example    : my $setsize = $array->setsize();
  Description: Getter, setter and lazy loader of setsize attribute for
               OligoArray objects. The setsize is the number of probes in a
               probeset for this array. This is likely to be 1 for non-Affy
               arrays.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub setsize {
	my $self = shift;
	$self->{'setsize'} = shift if @_;
	if ( !exists $self->{'setsize'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'setsize'};
}

=head2 type

  Arg [1]    : (optional) string - the type (currently either AFFY or OLIGO)
               for this array
  Example    : my $type = $array->type();
  Description: Getter, setter and lazy loader of type attribute for OligoArray
               objects. Currently the type can be either AFFY or OLIGO
  Returntype : string
  Exceptions : Throws if argument isn't a valid type (currently AFFY or OLIGO)
  Caller     : General
  Status     : Medium Risk

=cut

sub type {
	my $self = shift;
	my $type = shift;
	if ($type) {
		if ($VALID_TYPE{$type}) {
			$self->{'type'} = $type;
		} else {
			throw("$type is not a valid type for a Bio::EnsEMBL::OligoArray");
		}
	}
	if ( !exists $self->{'type'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'type'};
}

1;

