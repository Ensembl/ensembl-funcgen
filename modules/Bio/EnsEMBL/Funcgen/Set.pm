#
# Ensembl module for Bio::EnsEMBL::Funcgen::Set
#

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Funcgen::Set - A module to represent a base Set object.

=head1 SYNOPSIS

  use base qw( Bio::EnsEMBL::Funcgen::Set )

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
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );

use base qw( Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  MANDATORY ARGS:
  Arg [-NAME]          : String - name for this Set.
  Arg [-FEATURE_TYPE]  : Bio::EnsEMBL::Funcgen::FeatureType

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

  my ($name, $anal, $ftype, $ctype)
    = rearrange(['NAME', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE'], @_);

  #MANDATORY PARAMS
  throw('Need to specify a name')     if ! defined $name;
  assert_ref($ftype, 'Bio::EnsEMBL::Funcgen::FeatureType', 'Set FeatureType'); 
  assert_ref($anal, 'Bio::EnsEMBL::Analysis', 'Set Analysis'); 

  #OPTIONAL PARAMS
  if(defined $ctype){
    assert_ref($ctype, 'Bio::EnsEMBL::Funcgen::CellType', 'Set CellType'); 
  }

  #Define set_type automatically
  my @namespace = split/\:\:/, ref($self);
  ($self->{_set_type} = lc($namespace[$#namespace])) =~ s/set//;

  #Direct assignment as we have already validated
  $self->{name}         = $name;
  $self->{cell_type}    = $ctype;
  $self->{feature_type} = $ftype;
  $self->{analysis}     = $anal;

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

sub name { return shift->{name}; }


=head2 cell_type

  Example    : my $ctype_name = $set->cell_type->name;
  Description: Getter for the CellType for this Set.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub cell_type { return shift->{cell_type}; }


=head2 feature_type

  Example    : my $ftype_name = $set->feature_type->name;
  Description: Getter for the FeatureType of this Set.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_type { return shift->{feature_type}; }


=head2 analysis

  Example    : my $analysis_name = $set->analysis->logic_name;
  Description: Getter for the Analysis attribute of a Set.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub analysis {  return shift->{analysis}; }


=head2 set_type

  Example    : my $set_type = $set->set_type;
  Description: Getter for the set type attribute of this Set e.g. result, feature, input
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub set_type { return shift->{_set_type}; }



1;

