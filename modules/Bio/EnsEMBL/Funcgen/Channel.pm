#
# Ensembl module for Bio::EnsEMBL::Funcgen::Channel
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::Channel - A module to represent a single channel of an ExperimentalChip

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Channel;

my $array = Bio::EnsEMBL::Funcgen::Channel->new(
                   			        -EXPERIMENTAL_CHIP_ID => $ec_id,
			        	        -SAMPLE_ID            => $sample_id,
					        -TYPE                 => $type,
						-DYE                  => $dye,
						-DESCRIPTION          => $desc,
						);

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $chan_a = $db_adaptor->get_ChannelAdaptor();
my $chan = $chan_a->fetch_by_type_ExperimentalChip($type, $ExpChip);

=head1 DESCRIPTION

A Channel object represents a single channel on an ExperimentalChip. The data
are stored in the channel table, and associated expermental variables are
stored in the experimental_variable table.



=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Funcgen::Channel;


use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-EXPERIMENTAL_CHIP_ID]:  int - the experimental chip dbID



  Example    : my $array = Bio::EnsEMBL::Funcgen::Channel->new(
														       -EXPERIMENTAL_CHIP_ID => $ec_id,
 											                   -SAMPLE_ID            => $sample_id,
											                   -TYPE                 => $type,
											                   -DYE                  => $dye,
											                   -DESCRIPTION          => $desc,
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
	my ($ec_id, $sample_id, $type, $dye, $desc, $cell_line, $cell_line_id)
		= rearrange( ['EXPERIMENTAL_CHIP_ID', 'SAMPLE_ID', 'TYPE', 'DYE', 'DESCRIPTION', 'CELL_LINE', 'CELL_LINE_ID'], @_ );
	
	$self->sample_id($sample_id)          if defined $sample_id;
	$self->experimental_chip_id($ec_id)          if defined $ec_id;
	$self->type($type)    if defined $type;
	$self->dye($dye)    if defined $dye;
	$self->description($desc)   if defined $desc;
	$self->cell_line($cell_line) if defined $cell_line;
	$self->cell_line_id($cell_line_id) if defined $cell_line_id;

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

=head2 cell_line

  Arg [1]    : (optional) Bio::EnsemBL::Funcgen::CellLine
  Example    : my $cell_line_id = $chan->cell_line->dbID();
  Description: Getter, setter and lazy loader of the cell line attribute.
  Returntype : string
  Exceptions : Throws if arg not Bio::EnsEMBL::Funcgen::CellLine
  Caller     : General
  Status     : At Risk

=cut

sub cell_line {
  my $self = shift;

  if(@_){
    throw("Need to pass a valid Bio::EnsEMBL::Funcgen::CellLine object") if (! $_[0]->isa("Bio::EnsEMBL::Funcgen::CellLine"));
  }else{
    $self->{'cell_line'} = shift;
  }

  #this need to create an object from the id
  if ( ! exists $self->{'cell_line'} && $self->cell_line_id() && $self->adaptor() ) {
    $self->{'cell_line'} = $self->adaptor->db->get_CellLineAdaptor->fetch_by_dbID($self->cell_line_id());
  }

  return $self->{'cell_line'};
}

=head2 cell_line_id

  Arg [1]    : (optional) int - cell_line_id
  Example    : $chan->cell_line_id("1");
  Description: Getter, setter and lazy loader of the cell line id attribute.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub cell_line_id{
  my $self = shift;
  $self->{'cell_line_id'} = shift if @_;

  if ( ! exists $self->{'cell_line_id'} && $self->dbID() && $self->adaptor() ) {
    $self->adaptor->fetch_attributes($self);
  }

  return $self->{'cell_line_id'};
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
  
  if ( !exists $self->{'dye'} && $self->dbID() && $self->adaptor() ) {
    $self->adaptor->fetch_attributes($self);
  }
  return $self->{'dye'};
}


=head2 description

  Arg [1]    : (optional) string - the description of the Channel
  Example    : my $desc = $chan->description();
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

