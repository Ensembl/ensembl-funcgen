#
# Ensembl module for Bio::EnsEMBL::Funcgen::FeatureType
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::FeatureType - A module to represent a FeatureType.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::FeatureType;



=head1 DESCRIPTION


=head1 AUTHOR



This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::FeatureType;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new


=cut

sub new {
  my $class = shift;

  my $self = bless {},$class;
  
  
  my (
      $name,
      $desc,
     ) = rearrange([
		    'NAME', 'DESCRIPTION'
		   ], @_);
  
  
  if($name){
    $self->name($name);
  }else{
    throw("Must supply a FeatureType name\n");
  }

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
  Example    : my $desc = $probe->description();
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

1;

