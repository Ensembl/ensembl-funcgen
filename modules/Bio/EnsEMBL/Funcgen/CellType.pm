#
# Ensembl module for Bio::EnsEMBL::Funcgen::CellType
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::CellType - A module to represent a CellType.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::CellType;



=head1 DESCRIPTION

This is a simple class to represent information about a CellType.  This may represent an individual cell line or a more
generic tissue type.


=head1 AUTHOR

This module was written by Nathan Johnson.

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::CellType;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-name]          : string - name of CellType
  Arg [-display_label] : string - display label of CellType
  #Arg [-type]          : string - type of cell i.e.  ???? or Tissue, LINE?? enum?
  ?? xref to coriell ??
  Example              : my $ct = Bio::EnsEMBL::Funcgen::CellType->new(
                                                               -name  => "U2OS",
                                                               #-type  => "TISSUE",
                                                               -display_label => "Human Bone Osteosarcoma Epithelial Cells (U2OS)",
                                                                );
  Description: Constructor method for CellType class
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : Throws if name and type not defined.
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  
  
  my (
      $name,
      $dlabel,
     ) = rearrange([
		    'NAME', 'DISPLAY_LABEL',
		   ], @_);
  
  

  throw("Must supply a CellType name\n") if ! $name;

  $self->name($name);
  $self->display_label($dlabel) if $dlabel;

  return $self;
}



=head2 name

  Arg [1]    : string - name
  Example    : my $name = $ft->name();
  Description: Getter and setter of name attribute for CellType
               objects
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub name {
    my $self = shift;
    $self->{'name'} = shift if @_;
    return $self->{'name'};
}

=head2 display_label

  Arg [1]    : (optional) string - description
  Example    : my $display_label = $ct->display_label();
  Description: Getter and setter of display_label attribute for CellType
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
    my $self = shift;
    $self->{'display_label'} = shift if @_;
    return $self->{'display_label'};
}


=head2 class

  Arg [1]    : (optional) string - class
  Example    : $ft->class('HISTONE');
  Description: Getter and setter of description attribute for FeatureType
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

#sub class{
#  my $self = shift;
# $self->{'class'} = shift if @_;
#  return $self->{'class'};
#}
1;

