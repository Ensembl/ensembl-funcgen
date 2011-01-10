#
# Ensembl module for Bio::EnsEMBL::Funcgen::ArrayChip
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

Bio::EnsEMBL::Funcgen::ArrayChip - A simple module to represent the concept/template of 
a chip/slide within an array,  of which the physical manifestation is an ExperimentalChip.

=head1 SYNOPSIS

 use Bio::EnsEMBL::Funcgen::ArrayChip;

 my $ec = Bio::EnsEMBL::Funcgen::ArrayChip->new(
                    							 -ARRAY_ID  => $array->dbID(),
					                    		 -NAME      => $desc,
                                                 -DESIGN_ID => $design_id,
							                   );

#add more methods here?


=head1 DESCRIPTION

An ArrayChip object represent the concept of an array chip/slide withing a given array/chipset.
The data for ArrayChips is stored in the array_chip table.


=cut

use strict;
use warnings;


package Bio::EnsEMBL::Funcgen::ArrayChip;


use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-ARRAY_ID]  : int - the dbID of the parent array
  Arg [-ARRAY]     : Bio::EnsEMBL::Funcgen::Array
  Arg [-DESIGN_ID] : string - the unqiue deisng ID defined by the array vendor
  Arg [-NAME]      : string - the name of the array chip


  Example    : my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
							 -ARRAY_ID  => $array->dbID(),
							 -NAME      => $desc,
                                                         -DESIGN_ID => $design_id,
							 );								       );
  Description: Creates a new Bio::EnsEMBL::Funcgen::ArrayChip object.
  Returntype : Bio::EnsEMBL::Funcgen::ArrayChip
  Exceptions : None ? should throw if mandaotry params not set
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($array_id, $name, $design_id, $array)
    = rearrange( ['ARRAY_ID', 'NAME', 'DESIGN_ID', 'ARRAY'], @_ );

   
  #Remove array_id so we can remove checking below?

  throw("Must define a name($name) and design_id($design_id)") if(! $name || ! $design_id);


  #Make these mutually exclusive to avoid checking
  if($array_id && $array){
	throw('Must provide either -array or -array_id but not both');
  }

  if(defined $array){
	
	if(!(ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array'))){
	  throw('array paramter must be a valid Bio::EnsEMBL::Funcgen::Array');
	}
	
	$self->{'array'} = $array;
  }



  $self->array_id($array_id)  if defined $array_id;
  $self->name($name);
  $self->design_id($design_id);

  return $self;
}


=head2 array_id

  Arg [1]    : (optional) int - the parent array dbID
  Example    : my $array_id = $array_chip->array_id();
  Description: Getter, setter array_id attribute.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub array_id {
  my $self = shift;
  $self->{'array_id'} = shift if @_;
  return $self->{'array_id'};
}

=head2 name

  Arg [1]    : (optional) string - the array chip name
  Example    : my $ac_name = $array_chip->name();
  Description: Getter, setter for the name attribute
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if @_;
  return $self->{'name'};
}

=head2 design_id

  Arg [1]    : (optional) string - the array_chip unique design id as deinfed by the array vendor
  Example    : my $design_id = $array_chip->design_id();
  Description: Getter, setter for the design_id attribute
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub design_id {
  my $self = shift;
  $self->{'design_id'} = shift if @_; 
  return $self->{'design_id'};
}


=head2 get_Array

  Example    : my $array = $array_chip->get_array();
  Description: Getter for the array attribute
  Returntype : Bio::EnsEMBL::Funcgen::Array
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Array {
  my $self = shift;

  if(! defined $self->{'array'}){
    $self->{'array'} = $self->adaptor->db->get_ArrayAdaptor()->fetch_by_dbID($self->array_id());
  }

  return $self->{'array'};
}



1;

