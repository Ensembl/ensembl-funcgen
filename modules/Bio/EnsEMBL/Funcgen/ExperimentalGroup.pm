#
# Ensembl module for Bio::EnsEMBL::Funcgen::ExperimentalGroup
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

Bio::EnsEMBL::Funcgen::ExperimentalGroup - A module to represent 
an ExperimentalGroup. i.e. the authors of an experiment.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExperimentalGroup;



=head1 DESCRIPTION

This is a simple class to represent information about an ExperimentalGroup, 
containing a name and a more detailed description
This module is part of the Ensembl project: http://www.ensembl.org/

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::ExperimentalGroup;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-name]: string - name of ExperimentalGroup
  Arg [-location]: (optional) string - location of ExperimentalGroup
  Arg [-contact]: (optional) string - contact of ExperimentalGroup
  Arg [-url]: (optional) string - url containing information for the ExperimentalGroup
  Arg [-description]: (optional) string - descriptiom of ExperimentalGroup
  Arg [-project]: (optional) boolean - True if this is part of a large project (eg. ENCODE)
  Example    : my $group = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
                                                               -name  => "EBI",
                                                               -location  => "Hinxton",
                                                               -contact  => "dev@ensembl.org",
                                                               -url  => "http://www.ebi.ac.uk/",
                                                               -description => "European Bioinformatics Institute",
                                                               -is_project => 0,
                                                             );
  Description: Constructor method for ExperimentalGroup class
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalGroup
  Exceptions : Throws if name not defined 
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;

  my $obj_class = ref($caller) || $caller;
  my $self = $obj_class->SUPER::new(@_);
  
  my (
      $name,
      $location,
      $contact,
      $url,
      $desc,
      $is_project
     ) = rearrange([
					'NAME', 'LOCATION', 'CONTACT', 'URL', 'DESCRIPTION', 'IS_PROJECT'
				   ], @_);
  
  
  if($name){
    $self->name($name);
  }else{
    throw("Must supply a Group name\n");
  }
  $self->location($location) if $location;
  $self->contact($contact) if $contact;
  $self->url($url) if $url;
  $self->description($desc) if $desc;
  $self->is_project($is_project) if $is_project;

  return $self;
}



=head2 name

  Arg [1]    : string - name
  Example    : my $name = $ft->name();
  Description: Getter and setter of name attribute for ExperimentalGroup objects
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
  Example    : my $desc = $group->description();
  Description: Getter and setter of description attribute for ExperimentalGroup objects.
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

=head2 location

  Arg [1]    : (optional) string - location
  Example    : my $location = $group->location();
  Description: Getter and setter of location attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub location {
    my $self = shift;
    $self->{'location'} = shift if @_;
    return $self->{'location'};
}


=head2 contact

  Arg [1]    : (optional) string - contact
  Example    : my $contact = $group->contact();
  Description: Getter and setter of contact attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub contact {
    my $self = shift;
    $self->{'contact'} = shift if @_;
    return $self->{'contact'};
}


=head2 url

  Arg [1]    : (optional) string - url
  Example    : my $url = $group->url();
  Description: Getter and setter of url attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub url {
    my $self = shift;
    $self->{'url'} = shift if @_;
    return $self->{'url'};
}

=head2 is_project

  Arg [1]    : (optional) Boolean - is_project
  Example    : $group->is_project();
  Description: Getter and setter of is_project attribute for ExperimentalGroup objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : High Risk

=cut

sub is_project {
    my $self = shift;
    $self->{'is_project'} = shift if @_;
    return $self->{'is_project'};
}

1;

