#
# Ensembl module for Bio::EnsEMBL::Funcgen::Set
#

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Funcgen::Set - A module to represent a base Set object.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Funcgen::Set;

  @INC = qw (Bio::EnsEMBL::Funcgen::Set)

  sub new {
    my $caller = shift;
	
    my $class = ref($caller) || $caller;
	
    my $self = $class->SUPER::new(@_);

  
  }

=head1 DESCRIPTION

A base Set object which provides common methods available across Funcgen Set classes.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::FeatureSet
Bio::EnsEMBL::Funcgen::ResultSet
Bio::EnsEMBL::Funcgen::InputSet
Bio::EnsEMBL::Funcgen::Storable

=cut


package Bio::EnsEMBL::Funcgen::Set;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  MANDATORY ARGS:
  Arg [-NAME]          : String - name for this Set.
  Arg [-FEATURE_TYPE]  : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [-FEATURE_CLASS] : String - Class of feature e.g. result, annotated, 
                         regulatory, segmentation, external or dna_methylation.
  OPTIONAL ARGS:
  Arg [-CELL_TYPE]     : Bio::EnsEMBL::Funcgen::CellType
  Arg [-ANALYSIS]      : Bio::EnsEMBL::Analysis
  Arg [-DBID]          : Int
  Arg [-ADAPTOR]       : Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor e.g. Input|Result|FeatureSetAdaptor.

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor for Set objects.
  Returntype : Bio::EnsEMBL::Funcgen::Set
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

#Remove -type param this when fully implemented
#Removed -set_type param as this is auto generated from the namespace.
#Change set_type to mandatory and pass from ineritors?
#is_stored (dbID) check or leave to adaptor?

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;	
  my $self = $class->SUPER::new(@_);

  my ($name, $anal, $ftype, $ctype, $fclass, $type)
    = rearrange(['NAME', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE',
                 'FEATURE_CLASS', 'TYPE'], @_);
  
  #MANDATORY PARAMS
  throw('Need to specify a name')     if ! defined $name;
  $self->{feature_class} = $fclass || $type;
  throw('Need to specify a -feature_class') if ! defined $self->{feature_class};
  #feature_class validation should be done in inheritor constructors
  #against hash of valid enum field values

  if(! (ref($ftype) && $ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType')) ){
    throw('You must provide a valid Bio::EnsEMBL::Funcgen::FeatureType');
  }
  
  #OPTIONAL PARAMS
  if(defined $ctype && 
     ref($ctype) ne 'Bio::EnsEMBL::Funcgen::CellType'){
    throw('-CELL_TYPE param must be a valid Bio::EnsEMBL::Funcgen::CellType');
  }

  #Define set_type automatically
  my @namespace = split/\:\:/, ref($self);
  ($self->{_set_type} = lc($namespace[$#namespace])) =~ s/set//;	

  if(defined $anal){
    
    if(ref($anal) ne 'Bio::EnsEMBL::Analysis'){
      throw('-ANALYSIS argument must be a valid Bio::EnsEMBL::Analysis');
    }

    $self->{analysis} = $anal;
  }
  elsif($self->set_type ne 'input'){
    #Currently not mandatory for input_sets
    #Could move this to child Sets and just set analysis here
    #As with ftype
     throw('Must pass a valid -analysis parameter for a '.ref($self));
  }

  #Direct assignment as we have already validated
  $self->{name}         = $name;
  $self->{cell_type}    = $ctype;
  $self->{feature_type} = $ftype;
  
  return $self;
}


=head2 name

  Example    : my $set_name = $set->name;
  Description: Getter for the name of this Set.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return $_[0]->{name}; }


=head2 cell_type

  Example    : my $ctype_name = $set->cell_type->name;
  Description: Getter for the CellType for this Set.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub cell_type { return $_[0]->{cell_type}; }


=head2 feature_type

  Example    : my $ftype_name = $set->feature_type->name;
  Description: Getter for the FeatureType of this Set.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_type { return $_[0]->{feature_type}; }


=head2 feature_class

  Arg[0]     : String - feature class e.g. result, annotated, regulatory, external, dna_methylation or segmentation
  Example    : my $fclass = $set->feature_class;
  Description: Getter for the feature_type for this Set.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_class { return $_[0]->{feature_class}; }


=head2 feature_class_name

  Example    : my $fclass_adaptor_method = 'get_'.$set->feature_class.'Adaptor';
  Description: Getter for the full feature class name for this Set e.g. AnnotatedFeature, RegulatoryFeature
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

#Test for adaptor
#Can't set this in new as we won't pass the adaptor
#unless creating from _obj_from_sth

sub feature_class_name{
  my $self = shift;
  
  if(! defined $self->{feature_class_name} ){
    $self->{feature_class_name} = $self->adaptor->build_feature_class_name($self->feature_class);
  }
  
  return $self->{feature_class_name};
}



=head2 analysis

  Example    : my $analysis_name = $set->analysis->logic_name;
  Description: Getter for the Analysis attribute of a Set.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

#Not currently present in input_set
#implement in input_set instead of format

sub analysis {  return $_[0]->{analysis}; }


=head2 set_type

  Example    : my $set_type = $set->set_type;
  Description: Getter for the set type attribute of this Set e.g. result, feature, input
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub set_type { return $_[0]->{_set_type}; }




### DEPRECATED METHODS ###


=head2 type

  Example    : my $type = $set->type;
  Description: DEPRECATED Getter for the type for this Set.
               e.g. annotated, external, regulatory for FeatureSets
                    or 
                    array, sequencing for InputSets
               Currently not applicable to DataSets or ResultSets
  Exceptions : None
  Returntype : string 
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub type {
  my $self = shift;
  deprecate('Please use feature_class instead');
  return $self->feature_class(@_);
}



1;

