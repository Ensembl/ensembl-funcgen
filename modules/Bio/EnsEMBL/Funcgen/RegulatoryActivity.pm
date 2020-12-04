=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

  Bio::EnsEMBL::Funcgen::RegulatoryActivity

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
  );

  my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryFeature');
  my $regulatory_feature = $regulatory_feature_adaptor->fetch_by_stable_id('ENSR00000000011');

  my $regulatory_activity_list = $regulatory_activity_adaptor->fetch_all_by_RegulatoryFeature($regulatory_feature);

  print "The regulatory feature with stable id: "  . $regulatory_feature->stable_id . " has the following activities: \n";

  foreach my $current_regulatory_activity (@$regulatory_activity_list) {
    print "The activity in the epigenome "
      . $current_regulatory_activity->get_Epigenome->short_name
      . ' is: '
      . $current_regulatory_activity->activity
      . "\n";
  }

=head1 DESCRIPTION

=head1 SEE ALSO

=cut

package Bio::EnsEMBL::Funcgen::RegulatoryActivity;

use strict;
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_fetch
  _generic_set
);
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID                  => 'dbID',
    adaptor               => 'adaptor',
    adaptor               => 'adaptor',
    regulatory_feature_id => 'regulatory_feature_id',
    activity              => 'activity',
    epigenome_id          => 'epigenome_id',
  };
}

sub dbID                  { return shift->_generic_get_or_set('dbID', @_); }
sub regulatory_feature_id { return shift->_generic_get_or_set('regulatory_feature_id', @_); }
sub activity              { return shift->_generic_get_or_set('activity', @_); }
sub epigenome_id          { return shift->_generic_get_or_set('epigenome_id', @_); }

sub adaptor { 
  
  my $self = shift;

  my $adaptor = $self->_generic_get_or_set('adaptor', @_);
  
  if (defined $adaptor) {
    return $adaptor;
  }
  
  $adaptor = $self->get_RegulatoryFeature->adaptor;
  
  if (defined $adaptor) {
    return $adaptor;
  }
}

sub get_Epigenome {
  return shift->_generic_fetch('epigenome', 'get_EpigenomeAdaptor', 'epigenome_id');
}

sub set_Epigenome {
  my ($self, $obj) = @_;
  return $self->_generic_set('epigenome', 'Bio::EnsEMBL::Epigenome', $obj);
}

=head2 get_RegulatoryFeature

  Example       : say "Activity: "  .  $regulatory_activity->get_RegulatoryFeature;
  Description   : The regulatory feature this activity belongs to.
  Returns       : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Status        : At risk

=cut

sub get_RegulatoryFeature {
  return shift->_generic_fetch('regulatoryfeature', 'get_RegulatoryFeatureAdaptor', 'regulatory_feature_id');
}

sub set_RegulatoryFeature {
  my ($self, $obj) = @_;
  
  $self->_generic_set('regulatoryfeature', 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', $obj);
  use Scalar::Util qw( weaken );
  # Avoid circular reference and the memory leak that comes with it
  weaken($self->{regulatory_feature});
  return;

}

=head2 get_RegulatoryEvidence

  Example       : my $all_regulatory_evidence = $self->get_RegulatoryEvidence;
  Description   : Returns the regulatory evidence supporting this regulatory
                  activity. Regulatory evidence can be MotifFeatues
                  or AnnotatedFeatures.
  Returns       : ArrayRef[Bio::EnsEMBL::Funcgen::MotifFeature or Bio::EnsEMBL::Funcgen::AnnotatedFeature]
  Status        : At risk

=cut


sub get_RegulatoryEvidence {
  my ($self) = @_;

  if (defined $self->{'_regulatory_evidence'}) {
    return $self->{'_regulatory_evidence'}
  }
  my $regulatory_evidence_link = $self->get_RegulatoryEvidenceLink;

  if (! defined $regulatory_evidence_link) {
    return;
  }
  if (@$regulatory_evidence_link == 0) {
    return;
  }

  my $regulatory_feature = $self->get_RegulatoryFeature;
  my $regulatory_feature_slice = $regulatory_feature->slice;

  my @regulatory_evidence = map {
    $_->get_Evidence_on_Slice($regulatory_feature_slice)
  } @$regulatory_evidence_link;

  $self->{'_regulatory_evidence'} = \@regulatory_evidence;
  return $self->{'_regulatory_evidence'};
}


sub get_RegulatoryEvidenceLink {
  my ($self) = @_;

  if (defined $self->{'_regulatory_evidence_link'}) {
    return $self->{'_regulatory_evidence_link'}
  }

  $self->{'_regulatory_evidence_link'} = $self
    ->adaptor
    ->db
    ->get_RegulatoryEvidenceLinkAdaptor
    ->fetch_all_by_RegulatoryActivity($self);

    return $self->{'_regulatory_evidence_link'};
}

=head2 get_RegulatoryEvidence_by_type

  Arg [1]       : (mandatory) String, currently only 'motif'
  Example       : my $motif_evidence = $regulatory_activity->get_RegulatoryEvidence_by_type('motif');
  Description   : Returns links to regulatory evidence supporting this
                  regulatory activity..
  Returns       : ArrayRef[Bio::EnsEMBL::Funcgen::MotifFeature]
  Status        : At risk

=cut


sub get_RegulatoryEvidence_by_type {
  my ($self, $type) = @_;


  my $all_regulatory_evidence = $self->get_RegulatoryEvidence;
  my @evidence_by_type;

  if ($type eq 'motif') {
    @evidence_by_type = grep { $_->isa('Bio::EnsEMBL::Funcgen::MotifFeature')  } @$all_regulatory_evidence;
  }
  # if ($type eq 'annotated') {
  #   @evidence_by_type = grep { $_->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature')  } @$all_regulatory_evidence;
  # }
  return \@evidence_by_type;
}

=head2 feature_so_acc

  Example       : print $regulatory_activity->feature_so_acc;
  Description   : Returns the sequence ontology accession for this type of
                  regulatory feature.
  Returntype    : String
  Status        : At risk

=cut

sub feature_so_acc {
  my $self = shift;

  return $self->get_RegulatoryFeature->get_FeatureType->so_accession;
}

=head2 feature_so_term

  Example       : print $regulatory_activity->feature_so_term;
  Description   : Returns the sequence ontology term for this type of
                  regulatory feature.
  Returntype    : String
  Status        : At risk

=cut

sub feature_so_term {
  my $self = shift;

  return $self->get_RegulatoryFeature->get_FeatureType->so_term;
}

sub summary_as_hash {
  my ($self) = @_;

  my $regulatory_feature = $self->get_RegulatoryFeature;
  my $epigenome          = $self->get_Epigenome;
  my $feature_type       = $regulatory_feature->feature_type;

  return {
    regulatory_feature_stable_id => $regulatory_feature->stable_id,
    epigenome           => $epigenome->short_name,
    source              => $regulatory_feature->analysis->logic_name,
    bound_start         => $regulatory_feature->bound_seq_region_start,
    bound_end           => $regulatory_feature->bound_seq_region_end,
    start               => $regulatory_feature->seq_region_start,
    end                 => $regulatory_feature->seq_region_end,
    strand              => $regulatory_feature->strand,
    seq_region_name     => $regulatory_feature->seq_region_name,
    activity            => $self->activity,
    description         => $feature_type->name . ' ' . lc($self->activity) . ' in ' . $epigenome->short_name,
    feature_type        => $feature_type->name,
  };
}




1;
