#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSubset
#
# You may distribute this module under the same terms as Perl itself

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

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

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

