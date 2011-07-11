#
# Ensembl module for Bio::EnsEMBL::Funcgen::FeatureType
#


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

Bio::EnsEMBL::Funcgen::FeatureType - A module to represent a FeatureType. i.e. the target of an experiment.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::FeatureType;



=head1 DESCRIPTION

This is a simple class to represent information about a FeatureType, containing the name i.e Brno nomenclature or other controlled/validated name relevant to the class (HISTONE, PROMOTER etc), and description. This module is part of the Ensembl project: http://www.ensembl.org/

=cut

#To do
# add coding_transcript/gene methods.  Store as xrefs or custom feature_type_coding table? (miRanda etc)
# 

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::FeatureType;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-name]: string - name of FeatureType
  Arg [-class]: string - class of FeatureType
  Arg [-description]: string - descriptiom of FeatureType
  Example    : my $ft = Bio::EnsEMBL::Funcgen::FeatureType->new(
                                                               -name  => "H3K9Me",
                                                               -class => "HISTONE",
                                                               -description => "Generalised methylation of Histone 3 Lysine 9",
                                                                );
  Description: Constructor method for FeatureType class
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : Throws if name not defined ? and class
  Caller     : General
  Status     : Medium risk

=cut

sub new {
  my $caller = shift;

  my $obj_class = ref($caller) || $caller;
  my $self = $obj_class->SUPER::new(@_);
  
  my (
      $name,
      $desc,
      $class,
     ) = rearrange([
		    'NAME', 'DESCRIPTION', 'CLASS',
		   ], @_);
  
  
  if($name){
    $self->name($name);
  }else{
    throw("Must supply a FeatureType name\n");
  }


  #add test for class and enum? Validate names against Brno etc?
  $self->class($class) if $class;
  $self->description($desc) if $desc;

  return $self;
}



=head2 name

  Arg [1]    : string - name
  Example    : my $name = $ft->name();
  Description: Getter and setter of name attribute for FeatureType
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

=head2 description

  Arg [1]    : (optional) string - description
  Example    : my $desc = $ft->description();
  Description: Getter and setter of description attribute for FeatureType
               objects.
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

sub class{
  my $self = shift;
  $self->{'class'} = shift if @_;
  return $self->{'class'};
}


=head2 evidence_type_label

  Example    : my $track_label = $fsets[0]->feature_type->evidence_type_label.' MultiCell';
  Description: Getter for short evidence type label used in track label and field headers etc.
  Returntype : string
  Exceptions : None
  Caller     : Web code
  Status     : At risk

=cut

sub evidence_type_label{
  my $self = shift;

  #Could get undef key warn here
  #But only used in webcode so omit for speed

  return $Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor::regulatory_evidence_info{$self->class}->{label};
}

=head2 evidence_type_name

  Example    : my $long_name = $fsets[0]->feature_type->evidence_type_name;
  Description: Getter for evidence type name used in browser.
  Returntype : string
  Exceptions : None
  Caller     : Web code
  Status     : At risk

=cut

sub evidence_type_name{
  my $self = shift;

  #Could get undef key warn here
  #But only used in webcode so omit for speed

  return $Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor::regulatory_evidence_info{$self->class}->{name};
}




1;

