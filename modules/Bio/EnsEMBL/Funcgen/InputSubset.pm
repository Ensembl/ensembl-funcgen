#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSubset
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

Bio::EnsEMBL::Funcgen::InputSubset - A module to represent InputSubset object.


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::InputSubset;

my $input_subset = Bio::EnsEMBL::Funcgen::InputSubset->new
                    (
                     -NAME        => $name,
                     -INPUT_SET   => $iset,
                     -archive_id  => $archive_id,
                     -display_url => $display_url,
                     -replicate   => $iss_rep,
                     -is_control  => $is_control,
                    );

($input_subset) = @{$input_subset_adaptor->store($input_subset)};


my $control = ($input_set->is_control) ? 'control' : ''; 
  
print 'InputSubset ('.$input_subset->archive_id.
  ") is $control replicate ".$input_set->replicate.
  ":\t".$input_set->name."\nViewable here:\t".$input_set->display_url."\n";



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
                            (
                             -cell_type     => $cell_type,
                             -experiment    => $exp,
                             -feature_type  => $feature_type,
                             -archive_id    => $archive_id,
                             -display_url   => $display_url,
                             -is_control    => $is_control,
                             -name          => $name,
                             -replicate     => $iss_rep,
                            );


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

  my (
      $cell_type,
      $exp,
      $feature_type,
      $archive_id,
      $display_url,
      $is_control,
      $name,
      $rep,
      $analysis
      )
    = rearrange([
        'CELL_TYPE',
        'EXPERIMENT',
        'FEATURE_TYPE',
        'ARCHIVE_ID',
				'DISPLAY_URL',
        'IS_CONTROL',
        'NAME',
        'REPLICATE',
        'ANALYSIS',
        ], @_);

  throw('Must provide a name argument') if !defined $name;

  assert_ref($exp,          'Bio::EnsEMBL::Funcgen::Experiment');
  assert_ref($cell_type,    'Bio::EnsEMBL::Funcgen::CellType');
  assert_ref($feature_type, 'Bio::EnsEMBL::Funcgen::FeatureType');

  $self->{cell_type}    = $cell_type;
  $self->{experiment}   = $exp;
  $self->{feature_type} = $feature_type;
  $self->{archive_id}   = $archive_id;
  $self->{display_url}  = $display_url;
  $self->{is_control}   = $is_control;
  $self->{name}         = $name;
  $self->{replicate}    = $rep;

  return $self;
}


=head2 name

  Example    : my $name = $iss->name();
  Description: Getter for the name of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return shift->{name}; }


=head2 cell_type

  Example    : my $cell_type = $iss->cell_type;
  Description: Getter for the cell_type attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub cell_type { return shift->{cell_type}; }


=head2 experiment

  Example    : my $experiment = $iss->experiment;
  Description: Getter for the experiment attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub experiment { return shift->{experiment}; }


=head2 feature_type

  Example    : my $feature_type = $iss->feature_type;
  Description: Getter for the feature_type attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_type { return shift->{feature_type}; }


=head2 archive_id

  Example    : my $archive_id = $iss->archive_id;
  Description: Getter for the archive of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub archive_id { return shift->{archive_id}; }


=head2 display_url

  Example    : my $url = $iss->display_url;
  Description: Getter for the display_url of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_url{ return shift->{display_url}; }


=head2 replicate

  Example    : my $rep = $iss->replicate;
  Description: Getter for the replicate attribute of this InputSubset.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub replicate { return shift->{replicate}; }


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
                -experiment     => Bio::EnsEMBL::Funcgen::Experiment,

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

  my ($cell_type, $experiment, $feature_type) =
    rearrange(['CELL_TYPE', 'EXPERIMENT', 'FEATURE_TYPE'],
        %$params_hash);

  assert_ref($cell_type,    'Bio::EnsEMBL::Funcgen::CellType');
  assert_ref($experiment,   'Bio::EnsEMBL::Funcgen::Experiment');
  assert_ref($feature_type, 'Bio::EnsEMBL::Funcgen::FeatureType');

  $self->{cell_type}     = $cell_type;
  $self->{experiment}    = $experiment;
  $self->{feature_type}  = $feature_type;

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
               Defaults to: cell_type, experiment, feature_type
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

  $scl_methods ||= [qw(name archive_id display_url replicate is_control)];
  $obj_methods ||= [qw(cell_type experiment feature_type analysis)];

  return $self->SUPER::compare_to($obj, $shallow, $scl_methods, $obj_methods);
}


##### Deprecated ####

=head2 input_set

  Example    : my $input_set = $input_sset->input_set;
  Description: Getter for the input_set attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : None
  Caller     : General
  Status     : Deprecated, could return multiple InputSets

=cut

sub input_set {
  throw('v74, deprecated, please use InputSetAdaptor->fetch_all_by_InputSubsets');
}
1;

