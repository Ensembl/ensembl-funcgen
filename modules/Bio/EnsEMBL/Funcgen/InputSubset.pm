#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSubset
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

Bio::EnsEMBL::Funcgen::InputSubset - A module to represent InputSubset object.


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::InputSubset;

my $input_subset = Bio::EnsEMBL::Funcgen::InputSubset->new
                    (-epigenome            => $epigenome,
                     -experiment           => $exp,
                     -feature_type         => $feature_type,
                     -is_control           => $is_control,
                     -name                 => $name,
                     -biological_replicate => $iss_br,
                     -technical_replicate) => $iss_tr,;

($input_subset) = @{$input_subset_adaptor->store($input_subset)};

my $control = ($input_set->is_control) ? 'control' : ''; 
  
print 'InputSubset ('.$input_subset->name.
  ") is $control replicate ".$input_set->replicate."\n";

=head1 DESCRIPTION

An InputSubset object represents an individual distinct input within a given InputSet. This
normally translates to single file or replicate. There is no dedicated InputSubsetAdaptor,
store and fetch functionality is embedded within the InputSetAdaptor.

=cut

package Bio::EnsEMBL::Funcgen::InputSubset;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument   qw( rearrange );
use Bio::EnsEMBL::Utils::Exception  qw( throw deprecate );
use Bio::EnsEMBL::Utils::Scalar     qw( assert_ref );

use base qw( Bio::EnsEMBL::Funcgen::Set );

=head2 new

  Example    : my $iss = Bio::EnsEMBL::Funcgen::InputSubset->new
                            (-epigenome            => $epigenome,
                             -experiment           => $exp,
                             -feature_type         => $feature_type,
                             -is_control           => $is_control,
                             -name                 => $name,
                             -biological_replicate => $iss_br,
                             -technical_replicate  => $iss_tr,);

  Description: Constructor for InputSubset objects.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : Throws if no name defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($is_control, $name, $br, $tr) = rearrange
   (['IS_CONTROL', 'NAME', 'BIOLOGICAL_REPLICATE', 'TECHNICAL_REPLICATE'], @_);

  #FeatureType and Analysis validated in Set
  #Epigenome and Experiment validated in Set if defined  
 
  if(! defined $self->epigenome){
    throw('Mandatory parameter -epigenome is not defined');
  }

  if(! defined $self->experiment){
    throw('Mandatory parameter -experiment is not defined');
  }

  if(! defined $is_control){
    throw('Must defined an -is_control parameter'); 
    #is_control cannot be undef, as this will resolve to false
    #when storing
  }


  $self->{is_control}           = $is_control;
  $self->{name}                 = $name;
  $self->{biological_replicate} = $br;
  $self->{technical_replicate}  = $tr;

  #replicates are fine undef as they are not a boolean field

  return $self;
}


=head2 replicate

  Example    : my $rep = $iss->replicate;
  Description: Getter for the replicate attribute of this InputSubset.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Deprecated

=cut

sub replicate { 
  deprecate(
    'Bio::EnsEMBL::Funcgen::InputSubset::replicate has been deprecated and will be removed in Ensembl release 89 '
        . 'Please use Bio::EnsEMBL::Funcgen::InputSubset::technical_replicate instead.'
  );

  return shift->{technical_replicate};
}

=head2 biological_replicate

  Example    : my $br = $iss->biological_replicate;
  Description: Getter for the biological replicate attribute of this InputSubset.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub biological_replicate { return shift->{biological_replicate}; }

=head2 technical_replicate

  Example    : my $tr = $iss->technical_replicate;
  Description: Getter for the technical replicate attribute of this InputSubset.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub technical_replicate { return shift->{technical_replicate}; }


=head2 is_control

  Example    : if($iss->is_control){ # Do some control specific stuff here }
  Description: Getter for the is_control attribute of this InputSubset.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub is_control { return shift->{is_control}; }


=head2 reset_relational_attributes

  Arg[1]     : - Hashref containing the following parameters:
                -epigenome      => Bio::EnsEMBL::Funcgen::Epigenome,
                -experiment     => Bio::EnsEMBL::Funcgen::Experiment,
                -feature_type   => Bio::EnsEMBL::Funcgen::FeatureType,
                -analysis       => Bio::EnsEMBL::Funcgen::Analysis,

  Description: Resets all the relational attributes of a given InputSubset.
               Useful when creating a cloned object for migration beween DBs
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = @_;

  if(ref($params_hash) ne 'HASH'){
    throw('Must pass a HASHREF, not: ' .ref($params_hash));
  }

  my ($epigenome, $experiment, $feature_type, $analysis) =
    rearrange(['EPIGENOME', 'EXPERIMENT', 'FEATURE_TYPE', 'ANALYSIS'],
        %$params_hash);

  assert_ref($epigenome,    'Bio::EnsEMBL::Funcgen::Epigenome');
  assert_ref($experiment,   'Bio::EnsEMBL::Funcgen::Experiment');
  assert_ref($feature_type, 'Bio::EnsEMBL::Funcgen::FeatureType');
  assert_ref($analysis,     'Bio::EnsEMBL::Analysis');
  
  $self->{epigenome}     = $epigenome;
  $self->{experiment}    = $experiment;
  $self->{experiment_id} = $experiment->dbID;
  $self->{feature_type}  = $feature_type;
  $self->{analysis}      = $analysis;

  # Undef dbID and adaptor by default
  if(! $no_db_reset){
    $self->{dbID}    = undef;
    $self->{adaptor} = undef;
  }

  return;
}


=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::Storable (mandatory)
  Args[2]    : Boolean - Optional 'shallow' - no object methods compared
  Args[3]    : Arrayref - Optional list of InputSubset method names each
               returning a Scalar or an Array or Arrayref of Scalars.
               Defaults to: name archive_id display_url replicate is_control
  Args[4]    : Arrayref - Optional list of InputSubset method names each
               returning a Storable or an Array or Arrayref of Storables.
               Defaults to: epigenome, experiment, feature_type
  Example    : my %shallow_diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this InputSubset to another based on the defined scalar
               and storable methods.
  Returntype : Hashref of key attribute/method name keys and values which differ.
               Keys will always be the method which has been compared.
               Values can either be a error string, a hashref of diffs from a
               nested object, or an arrayref of error strings or hashrefs where
               a particular method returns more than one object.
  Exceptions : None
  Caller     : Import/migration pipeline
  Status     : At Risk

=cut


sub compare_to {
  my ($self, $obj, $shallow, $scl_methods, $obj_methods) = @_;

  $scl_methods ||= [qw(name replicate is_control)];
  $obj_methods ||= [qw(epigenome experiment feature_type analysis)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods, $obj_methods);
}


1;

