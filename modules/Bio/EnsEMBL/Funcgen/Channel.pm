#
# Ensembl module for Bio::EnsEMBL::Funcgen::Channel
#
# You may distribute this module under the same terms as Perl itself


=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Funcgen::Channel - A module to represent a single channel of an ExperimentalChip

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Channel;

my $array = Bio::EnsEMBL::Funcgen::Channel->new(
                   			        -EXPERIMENTAL_CHIP_ID => $ec_id,
			        	        -SAMPLE_ID            => $sample_id,
					        -TYPE                 => $type,
						-DYE                  => $dye,
					       );

#-replace TYPE with DENOMINATOR?

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $chan_a = $db_adaptor->get_ChannelAdaptor();
my $chan = $chan_a->fetch_by_type_ExperimentalChip($type, $ExpChip);

=head1 DESCRIPTION

A Channel object represents a single channel on an ExperimentalChip. The data
are stored in the channel table, and associated expermental variables are
stored in the experimental_variable table.

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Funcgen::Channel;


use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-EXPERIMENTAL_CHIP_ID]:  int - the experimental chip dbID



  Example    : my $array = Bio::EnsEMBL::Funcgen::Channel->new(
	                   			       -EXPERIMENTAL_CHIP_ID => $ec_id,
                                   	               -SAMPLE_ID            => $sample_id,
					               -TYPE                 => $type,
					               -DYE                  => $dye,
					                     );
  Description: Creates a new Bio::EnsEMBL::Funcgen::Channel object.
  Returntype : Bio::EnsEMBL::Funcgen::Channel
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);

	#can we lc these?
	my ($ec_id, $sample_id, $type, $dye)
		= rearrange( ['EXPERIMENTAL_CHIP_ID', 'SAMPLE_ID', 'TYPE', 'DYE'], @_ );
	
	$self->sample_id($sample_id)          if defined $sample_id;
	$self->experimental_chip_id($ec_id)          if defined $ec_id;
	$self->type($type)    if defined $type;
	$self->dye($dye)    if defined $dye;


	return $self;
}

=head2 sample_id

  Arg [1]    : (optional) string - the sample id for this Channel
  Example    : my $sample_id = $chan->sample_id();
  Description: Getter, setter and lazy loader of sample_id attribute.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub sample_id {
	my $self = shift;
	$self->{'sample_id'} = shift if @_;

	if ( ! exists $self->{'sample_id'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}

	return $self->{'sample_id'};
}



=head2 experimental_chip_id

  Arg [1]    : (optional) int - the experimenta chip dbID
  Example    : my $ec_id = $chan->experimental_chip_id();
  Description: Getter, setter and lazy loader of experimental_chip_id attribute
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub experimental_chip_id {
	my $self = shift;

	$self->{'experimental_chip_id'} = shift if @_;
	if ( !exists $self->{'experimental_chip_id'} && $self->dbID() && $self->adaptor() ) {
		$self->adaptor->fetch_attributes($self);
	}
	return $self->{'experimental_chip_id'};
}

=head2 type

  Arg [1]    : (optional) string - the channel type e.g. EXPERIMENTAL or CONTROL
  Example    : my $type = $chan->type();
  Description: Getter, setter and lazy loader of type attribute
  Returntype : string
  Exceptions : 
  Caller     : General
  Status     : Medium Risk

=cut

sub type {
  my $self = shift;

  $self->{'type'} = shift if @_;

  #warn "we need to control EXPERIMENTAL OR CONTROL here or enum on DB";
  
  if ( !exists $self->{'type'} && $self->dbID() && $self->adaptor() ) {
    $self->adaptor->fetch_attributes($self);
  }
  return $self->{'type'};
}

=head2 dye

  Arg [1]    : (optional) string - the channel type e.g. EXPERIMENTAL or CONTROL
  Example    : my $dye = $chan->dye();
  Description: Getter, setter and lazy loader of dye attribute
  Returntype : string
  Exceptions : 
  Caller     : General
  Status     : Medium Risk

=cut

sub dye {
  my $self = shift;

  $self->{'dye'} = shift if @_;
  
  #if ( !exists $self->{'dye'} && $self->dbID() && $self->adaptor() ) {
  #  $self->adaptor->fetch_attributes($self);
  #}
  return $self->{'dye'};
}


1;

