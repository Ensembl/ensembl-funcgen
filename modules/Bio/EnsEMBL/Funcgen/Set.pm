#
# Ensembl module for Bio::EnsEMBL::Funcgen::Set
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

A base Set object which provides common methods available across Funcgen Set classes. (Apart from DataSet).

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::FeatureSet
Bio::EnsEMBL::Funcgen::ResultSet
Bio::EnsEMBL::Funcgen::InputSubSet
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
  Arg [-ANALYSIS]      : Bio::EnsEMBL::Analysis

  OPTIONAL ARGS:
  Arg [-EPIGENOME]     : Bio::EnsEMBL::Funcgen::Epigenome
  Arg [-DBID]          : Int
  Arg [-ADAPTOR]       : Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor e.g. Result|Feature|InputSubSetAdaptor.

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

  my ( $name, $anal, $ftype, $epigenome, $exp, $exp_id )
    = rearrange( [ 'NAME',         'ANALYSIS',
                   'FEATURE_TYPE', 'EPIGENOME',
                   'EXPERIMENT',   'EXPERIMENT_ID' ],
                 @_ );

  #MANDATORY PARAMS
  throw('Need to specify a name')     if ! defined $name;

  use Carp;
  
  if (! defined $ftype) {
    confess("The feature type was undefined!");
  }
  
  if (! $ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType')) {
    confess("The feature type passed was a " . (ref $ftype));
  }
  
  assert_ref($anal, 'Bio::EnsEMBL::Analysis', 'Set Analysis');

  if($exp){
    assert_ref($exp, 'Bio::EnsEMBL::Funcgen::Experiment', 'Set Experiment');
    #Just ignoring/overwriting $exp_id here
    $exp_id = $exp->dbID;
  }

  #OPTIONAL PARAMS
  if(defined $epigenome){
    assert_ref($epigenome, 'Bio::EnsEMBL::Funcgen::Epigenome', 'Set Epigenome'); 
  }

  #Define set_type automatically
  my @namespace = split/\:\:/, ref($self);
  ($self->{_set_type} = lc($namespace[$#namespace])) =~ s/set//;

  #Direct assignment as we have already validated
  $self->{name}          = $name;
  $self->{epigenome}     = $epigenome;
  $self->{feature_type}  = $ftype;
  $self->{analysis}      = $anal;

  #Only set this if defined as we use exists as test
  #test for undef experiment
  $self->{experiment}    = $exp if $exp;

  $self->{experiment_id} = $exp_id;
  #We need the exp_id to enable quick source_label fetch
  #without having to create the Experiment

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


=head2 epigenome

  Example    : my $epigenome = $set->epigenome->name;
  Description: Getter for the Epigenome for this Set.
  Returntype : Bio::EnsEMBL::Funcgen::Epigenome
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub epigenome {
    return shift->{epigenome};
}


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



=head2 experiment_id

  Example    : my $exp_id = $set->experiment_id;
  Description: Getter for the experiment_id attribute of this Set 
  Returntype : Int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut


sub experiment_id { return shift->{experiment_id}; }



=head2 experiment

  Example    : my $exp = $set->experiment;
  Description: Getter for the Experiment of this Set.
               Returns undef if there is no experiment_id associatied with this Set
  Returntype : Bio::EnsEMBL::Funcgen::Experiment or undef
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub experiment{ 
  my $self = shift;
  
  if ((! exists $self->{experiment}) && 
      $self->experiment_id){
    $self->{experiment} = $self->adaptor->db->get_ExperimentAdaptor->fetch_by_dbID($self->experiment_id);
  }
     
  return $self->{experiment};
}

1;

