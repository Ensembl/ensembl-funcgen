#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSubset
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

Bio::EnsEMBL::InputSet - A module to represent InputSubset object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExperimetnalSubset;

my $data_set = Bio::EnsEMBL::Funcgen::InputSubset->new(
	                                                         -DBID             => $dbID,
							 					             -ADAPTOR          => $self,
                                                             -NAME             => $name,
                                                             -EXPERIMENTAL_SET => $eset,
                                                             );



=head1 DESCRIPTION

An InputSubset object is a very simple skeleton class to enable storage of associated subset states. As such there
are only very simple accessor methods for basic information, and there is no namesake adaptor, rather is is handled by the 
InputSetAdaptor.


=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::InputSubset;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Example    : my $eset = Bio::EnsEMBL::Funcgen::InputSubset->new(
                                                                        -DBID            => $dbID,
							 					                        -ADAPTOR         => $self,
                                                                        -NAME            => $name,
                                                                        -EXPERIMENTAL_SET => $eset,
                                                                        );


  Description: Constructor for InputSubset objects.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : Throws if no name defined
               Throws if CellType or FeatureType are not valid or stored
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  #do we need to add $fg_ids to this?  Currently maintaining one feature_group focus.(combi exps?)
  my ($name, $eset)
    = rearrange(['NAME', 'INPUT_SET'], @_);
  
  
  throw('Must provide a name argument') if ! defined $name;

  if(!(ref($eset) && 
	   $eset->isa('Bio::EnsEMBL::Funcgen::InputSet')
	   && $eset->dbID())){
	throw('Must provide a valid stored input_set argument');
  }
  

  $self->{'name'} = $name;
  $self->{'input_set'} = $eset;

  return $self;
}


=head2 name

  Example    : my $name = $exp_sset->name();
  Description: Getter for the name of this InputSubset.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub name {
  my $self = shift;
  return $self->{'name'};
}

=head2 input_set

  Example    : my $eset = $exp_sset->input_set();
  Description: Getter for the input_set attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub input_set {
  my $self = shift;
  return $self->{'input_set'};
}



1;

