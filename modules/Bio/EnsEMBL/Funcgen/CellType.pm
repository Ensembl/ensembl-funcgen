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

This is a simple class to represent information about a CellType.  
This may represent harvested cells, a cell line or a more generic tissue type.


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

  Arg [1]    : String - name of CellType
  Arg [2]    : String - display label of CellType
  Arg [3]    : String - description of CellType
  Arg [4]    : String - gender e.g. male, female or NULL
  Arg [5]    : String - Experimental Factor Ontology ID e.g. EFO_0002869
 
  Example              : my $ct = Bio::EnsEMBL::Funcgen::CellType->new
                                    (
                                     -name          => "U2OS",
                                     -display_label => "U20S",
                                     -description   => "Human Bone Osteosarcoma Epithelial Cells",
                                     -gender        => 'female',
                                     -efo_id        => 'EFO_0002869',
                                    );

  Description: Constructor method for CellType class
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : Throws if name not defined.
  Caller     : General
  Status     : Stable

=cut

#-type/class => "TISSUE", enum? Mandatory.
#remove display label?

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  
  
  my (
      $name,
      $dlabel,
	  $desc,
	  $gender,
	  $efo_id
     ) = rearrange([
		    'NAME', 'DISPLAY_LABEL', 'DESCRIPTION','GENDER', 'EFO_ID'
		   ], @_);
  
  

  throw("Must supply a CellType name") if ! $name;

  if(defined $gender){
	#enum will not force this so validate here
	throw("Gender must be either male or female") if ! grep(/^$gender$/, ('male', 'female'));
	$self->gender($gender);
  }

  $self->{'name'} = $name; #Set directly as mandatory, to enable getter only method
  $self->display_label($dlabel) if defined $dlabel;
  $self->description($desc)     if defined $desc;
  $self->efo_id($efo_id)        if defined $efo_id;

  return $self;
}



=head2 name

  Arg [1]    : String - name
  Example    : my $name = $ct->name();
  Description: Getter  of name attribute for CellType
               objects
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub name {
    my $self = shift;
    return $self->{'name'};
}

=head2 gender

  Arg [1]    : String (optional) - gender e.g. male or female
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

  Arg [1]    : String (optional) - description
  Example    : my $desc = $ct->description();
  Description: Getter and setter of description attribute for CellType
               objects
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

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
    return $self->{'display_label'};
}


=head2 efo_id

  Arg [1]    : (optional) String - Experimental Factor Ontology ID
  Example    : $ft->efo_id('EFO_0001196');
  Description: Getter and setter of the Experimental Factor Ontology ID
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub efo_id{
  my $self = shift;
  $self->{'efo_id'} = shift if @_;
  return $self->{'efo_id'};
}
1;

