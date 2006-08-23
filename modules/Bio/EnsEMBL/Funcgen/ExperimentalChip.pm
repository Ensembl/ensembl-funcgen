#
# Ensembl module for Bio::EnsEMBL::Funcgen::ExperimentalChip
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::Experimental - A module to represent a physical unique experimental chip.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExperimentalChip;

my $array = Bio::EnsEMBL::Funcgen::ExperimentalChip->new(
														 -dbID           => $ec_id,
														 -unique_id => $c_uid,
														 -experiment_id  => $exp_id,
	                                                     -array_chip_id  => $ac_id,
														 -description    => $desc,
														 );

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $ec_a = $db_adaptor->get_ExperimentalChipAdaptor();
my $echip = $ec_a->fetch_by_unique_and_experimental_id($c_uid, $exp_dbid)

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


package Bio::EnsEMBL::Funcgen::ExperimentalChip;


use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-unique_id]: int - the unique id of this individual experimental chip
  Arg [-experiment_id]:  int - the experiment dbID



  Example    : my $array = Bio::EnsEMBL::Funcgen::ExperimentalChip->new(
														 -dbID           => $ec_id,
														 -unique_id => $c_uid,
														 -experiment_id  => $exp_id,
                                                       	 -array_chip_id  => $ac_id,
														 -description    => $desc,
														 );
								 );
  Description: Creates a new Bio::EnsEMBL::Funcgen::ExperimentalChip object.
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);

	#can we lc these?
	my ($c_uid, $exp_dbid, $ac_id, $desc)
		= rearrange( ['UNIQUE_ID', 'EXPERIMENT_ID', 'ARRAY_CHIP_ID', 'DESCRIPTION'], @_ );
	
	$self->unique_id($c_uid)          if defined $c_uid;
	$self->experiment_id($exp_dbid)          if defined $exp_dbid;
	$self->array_chip_id($ac_id)    if defined $ac_id;
	$self->description($desc)   if defined $desc;


	return $self;
}

=head2 get_Channels

  Args       : None
  Example    : my $channels = $array->get_Channels();
  Description: Returns all channels on a ExperimentalChip. Needs a database connection.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_Channels {
	my $self = shift;
	
	if ( $self->dbID() && $self->adaptor() ) {
		foreach my $channel (@{$self->adaptor->db->get_ChannelAdaptor->fetch_all_by_ExperimentalChip($self)}){
			$self->{'channels'}{$channel->dbID()} = $channel;
		}	
	} else {
		warning('Need database connection to retrieve Channels');
	}

	return [values %{$self->{'channels'}}];
}




=head2 unique_id

  Arg [1]    : (optional) int - the unique chip id for this ExperimentalChip
  Example    : my $c_uid = $array->unique_id();
  Description: Getter, setter and lazy loader of unique_id attribute.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub unique_id {
	my $self = shift;
	$self->{'unique_id'} = shift if @_;

	if ( ! exists $self->{'unique_id'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}

	return $self->{'unique_id'};
}

=head2 experiment_id

  Arg [1]    : (optional) int - the experiment dbID
  Example    : my $exp_id = $array->experiment_id();
  Description: Getter, setter and lazy loader of experiment_id attribute
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub experiment_id {
	my $self = shift;

	$self->{'experiment_id'} = shift if @_;
	if ( !exists $self->{'experiment_id'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'experiment_id'};
}

=head2 array_chip_id

  Arg [1]    : (optional) int - the array_chip dbID
  Example    : my $ac_id = $ec->array_chip_id();
  Description: Getter, setter and lazy loader of array_chip_id attribute
  Returntype : int
  Exceptions : 
  Caller     : General
  Status     : Medium Risk

=cut

sub array_chip_id {
  my $self = shift;

  $self->{'array_chip_id'} = shift if @_; 

  if ( !exists $self->{'array_chip_id'} && $self->dbID() && $self->adaptor() ) {
    $self->adaptor->fetch_attributes($self);
  }
  return $self->{'array_chip_id'};
}


=head2 description

  Arg [1]    : (optional) string - the description of the ExperimentalChip
  Example    : my $desc = $ec->description();
  Description: Getter, setter and lazy loader of description attribute
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

