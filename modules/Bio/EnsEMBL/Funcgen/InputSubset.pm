#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSubset
#


=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

Bio::EnsEMBL::InputSubset - A module to represent InputSubset object.


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::InputSubset;

my $input_subset = Bio::EnsEMBL::Funcgen::InputSubset->new
                    (
                     -DBID        => $dbID,
                     -ADAPTOR     => $self,
                     -NAME        => $name,
                     -INPUT_SET   => $iset,
                     -archive_id  => $archive_id,
                     -display_url => $display_url,
                     -replicate   => $iss_rep,
                     -is_control  => $is_control,
                    );



=head1 DESCRIPTION

An InputSubset object represents an individual distinct input within a given InputSet. This
normally translates to single file or replicate. There is no dedicated InputSubsetAdaptor,
store and fetch functionality is embedded within the InputSetAdaptor.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::InputSubset;

use Bio::EnsEMBL::Utils::Argument   qw ( rearrange );
use Bio::EnsEMBL::Utils::Exception  qw ( throw deprecate );
use Bio::EnsEMBL::Utils::Scalar     qw ( assert_ref );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


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
      $exp,
      $archive_id,
      $display_url,
      $is_control,
      $name,
      $rep,
      )
    = rearrange([
        'EXPERIMENT',
        'ARCHIVE_ID',
				'DISPLAY_URL',
        'IS_CONTROL',
        'NAME',
        'REPLICATE',
        ], @_);

  throw('Must provide a name argument') if(!defined $name);

  assert_ref($exp, 'Bio::EnsEMBL::Funcgen::Experiment');
  #can't do is_stored_and_valid here as we don't adaptor
    if(!$exp->dbID){
    throw('Must provide a valid stored Bio::EnsEMBL::Funcgen::Experiment');
  }
# set in Set
  throw ('Must provide a  FeatureType') if(! defined $self->feature_type);
  throw ('Must provide a  CellType')    if(! defined $self->cell_type);

  $self->{experiment}  = $exp;
  $self->{archive_id}  = $archive_id;
  $self->{display_url} = $display_url;
  $self->{is_control}  = $is_control;
  $self->{name}        = $name;
  $self->{replicate}   = $rep;

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

sub name { return $_[0]->{name}; }


=head2 cell_type

  Example    : my $cell_type = $iss->cell_type;
  Description: Getter for the cell_type attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub cell_type {    return $_[0]->{cell_type}; }

=head2 experiment

  Example    : my $experiment = $iss->experiment;
  Description: Getter for the experiment attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub experiment {    return $_[0]->{experiment}; }

=head2 feature_type

  Example    : my $feature_type = $iss->feature_type;
  Description: Getter for the feature_type attribute of this InputSubset.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_type {    return $_[0]->{feature_type}; }


=head2 archive_id

  Example    : my $archive_id = $iss->archive_id;
  Description: Getter for the archive of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub archive_id { return $_[0]->{archive_id}; }


=head2 display_url

  Example    : my $url = $iss->display_url;
  Description: Getter for the display_url of this InputSubset.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_url{ return $_[0]->{display_url}; }


=head2 replicate

  Example    : my $rep = $iss->replicate;
  Description: Getter for the replicate attribute of this InputSubset.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

  #is_stored (in corresponding db) checks will be done in store method
sub replicate { return $_[0]->{replicate}; }


=head2 is_control

  Example    : if($iss->is_control){ # Do some control specific stuff here }
  Description: Getter for the is_control attribute of this InputSubset.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub is_control { return $_[0]->{is_control}; }

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
  $obj_methods ||= [qw(cell_type experiment feature_type)];

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

