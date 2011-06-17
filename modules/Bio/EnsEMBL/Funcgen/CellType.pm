#
# Ensembl module for Bio::EnsEMBL::Funcgen::CellType
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

Bio::EnsEMBL::Funcgen::CellType - A module to represent a CellType.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::CellType;



=head1 DESCRIPTION

This is a simple class to represent information about a CellType.  This may represent an individual cell line or a more
generic tissue type.


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
                                                               -name          => "U2OS",
                                                               -display_label => "",#?
                                                               -description   => "Human Bone Osteosarcoma Epithelial Cells",
                                                               #-type/class => "TISSUE", enum?
                                                               #xref/coriell id?
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
	  $desc,
	  $gender,
     ) = rearrange([
		    'NAME', 'DISPLAY_LABEL', 'DESCRIPTION','GENDER'
		   ], @_);
  
  

  throw("Must supply a CellType name") if ! $name;

  if(defined $gender){
	throw("Gender must be either male or female") if ! grep(/^$gender$/, ('male', 'female'));
	$self->gender($gender);
  }

  $self->name($name);
  $self->display_label($dlabel) if $dlabel;
  $self->description($desc) if $desc;

  return $self;
}



=head2 name

  Arg [1]    : string - name
  Example    : my $name = $ct->name();
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

=head2 gender

  Arg [1]    : string - gender e.g. male or female
  Example    : my $gender = $ct->gender();
  Description: Getter for the gender attribute for 
               CellType objects
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub gender {
    my $self = shift;
    $self->{'gender'} = shift if @_;
    return $self->{'gender'};
}

=head2 description

  Arg [1]    : string - description
  Example    : my $desc = $ct->description();
  Description: Getter and setter of description attribute for CellType
               objects
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if @_;
    return $self->{'description'};
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


	#This may result in very long dispay_labels
	#if(! defined $self->{'display_label'}){
	#	$self->{'display_label'} = (defined $self->{'desciption'}) ? $self->description()." (".$self->name().")" : $self->name();
	#}

    return $self->{'display_label'};
}


#=head2 class
#
#  Arg [1]    : (optional) string - class
#  Example    : $ft->class('HISTONE');
#  Description: Getter and setter of description attribute for FeatureType
#               objects.
#  Returntype : string
#  Exceptions : None
#  Caller     : General
#  Status     : Low Risk
#
#=cut

#sub class{
#  my $self = shift;
# $self->{'class'} = shift if @_;
#  return $self->{'class'};
#}
1;

