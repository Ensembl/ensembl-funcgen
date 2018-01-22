=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
      . $current_regulatory_activity->get_Epigenome->display_label 
      . ' is: ' 
      . $current_regulatory_activity->activity 
      . "\n";
  }

=head1 DESCRIPTION
=head1 SEE ALSO
=cut

package Bio::EnsEMBL::Funcgen::RegulatoryActivity;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self =  bless({}, $class);
  return $self;
}

sub db {
  my $self = shift;
  my $db   = shift;

  if ($db) {
    $self->{'_db'}  = $db;
  }
  return $self->{'_db'};
}

sub dbID {
  my $self = shift;
  my $dbID = shift;

  if ($dbID) {
    $self->{'_dbID'}  = $dbID;
  }
  return $self->{'_dbID'};
}

sub _epigenome_id {
  my $self = shift;
  my $epigenome_id = shift;

  if(defined $epigenome_id) {
    $self->{'_epigenome_id'}  = $epigenome_id;
  }
  return $self->{'_epigenome_id'};
}

=head2 get_Epigenome

  Description   : Deprecated, use get_Epigenome instead.

=cut

sub epigenome {
  my $self = shift;
  return $self->get_Epigenome;
}

=head2 get_Epigenome

  Example       : $regulatory_feature_summary = $regulatory_feature->get_Epigenome;
  Description   : Retrieves a textual summary of this RegulatoryFeature.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub get_Epigenome {
  my $self = shift;

  if(! defined $self->{'_epigenome'}) {
    my $ea = $self->db->{'_epigenome_adaptor'} ||= $self->db->get_EpigenomeAdaptor();
    $self->{'_epigenome'} = $ea->fetch_by_dbID($self->_epigenome_id);
  }
  return $self->{'_epigenome'};
}

sub _is_multicell {
  my $self = shift;
  my $is_multicell = shift;

  if(defined $is_multicell) {
    $self->{'_is_multicell'}  = $is_multicell;
  }
  return $self->{'_is_multicell'};
}

=head2 activity

  Example       : say "Activity: "  .  $regulatory_activity->activity;
  Description   : The activity this object represents as a string. 
  Returns       : String, one of 'INACTIVE','REPRESSED','POISED','ACTIVE' or 'NA'.
  Status        : At risk

=cut

sub activity {
  my $self = shift;
  my $activity = shift;

  if(defined $activity) {
    $self->{'_activity'}  = $activity;
  }
  return $self->{'_activity'};
}

sub _regulatory_feature_id {
  my $self = shift;
  my $regulatory_feature_id = shift;

  if(defined $regulatory_feature_id) {
    $self->{'_regulatory_feature_id'} = $regulatory_feature_id;
  }
  return $self->{'_regulatory_feature_id'};
}

=head2 get_RegulatoryFeature

  Example       : say "Activity: "  .  $regulatory_activity->get_RegulatoryFeature;
  Description   : The regulatory feature this activity belongs to.
  Returns       : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Status        : At risk

=cut

sub get_RegulatoryFeature {
  my $self = shift;

  if (defined $self->{'_regulatory_feature'}) {
    return $self->{'_regulatory_feature'}
  }
  if (defined $self->_regulatory_feature_id) {
  
    $self->{'_regulatory_feature'} = $self
      ->db
      ->get_RegulatoryFeatureAdaptor
      ->fetch_by_dbID(
	$self->_regulatory_feature_id
      );
      return $self->{'_regulatory_feature'};
  }
#   throw("Can't get regulatory feature. No regulatory feature object has been set nor a regulatory feature id by which to fetch one.");
  return undef;
}

=head2 get_RegulatoryEvidenceLink

  Example       : my $regulatory_evidence_link = $regulatory_activity->get_RegulatoryEvidenceLink;
  Description   : Returns links to regulatory evidence supporting this
                  regulatory activity.
  Returns       : ArrayRef[Bio::EnsEMBL::Funcgen::RegulatoryEvidenceLink]
  Status        : At risk

=cut

sub get_RegulatoryEvidenceLink {
  my $self = shift;

  if (defined $self->{'_regulatory_evidence_link'}) {
    return $self->{'_regulatory_evidence_link'}
  }

  $self->{'_regulatory_evidence_link'} = $self
    ->db
    ->get_RegulatoryEvidenceLinkAdaptor
    ->fetch_all_by_RegulatoryActivity($self);

    return $self->{'_regulatory_evidence_link'};
}

=head2 get_RegulatoryEvidence_by_type

  Arg [1]       : (mandatory) String, 'motif' or 'annotated'

  Example       : my $motif_evidence = $regulatory_activity->get_RegulatoryEvidence_by_type('motif');

  Description   : Returns links to regulatory evidence supporting this
                  regulatory activity.

                  There are two types of regulatory evidence: MotifFeatues 
                  and AnnotatedFeatures.

                  Depending on which type is needed, the first argument has 
                  to be either 'motif' or 'annotated'.

  Returns       : 
                  If Arg [1] is 'motif':     ArrayRef[Bio::EnsEMBL::Funcgen::MotifFeature]
                  If Arg [1] is 'annotated': ArrayRef[Bio::EnsEMBL::Funcgen::AnnotatedFeature]

  Status        : At risk

=cut

sub get_RegulatoryEvidence_by_type {
  my $self = shift;
  
  # 'motif' or 'annotated'
  my $type = shift;
  
  my $all_regulatory_evidence = $self->get_RegulatoryEvidence;
  my @evidence_by_type;

  if ($type eq 'motif') {
    @evidence_by_type = grep { $_->isa('Bio::EnsEMBL::Funcgen::MotifFeature')  } @$all_regulatory_evidence;
  }
  if ($type eq 'annotated') {
    @evidence_by_type = grep { $_->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature')  } @$all_regulatory_evidence;
  }
  return \@evidence_by_type;
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
  my $self = shift;

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

sub set_RegulatoryFeature {
  my $self = shift;
  my $regulatory_feature = shift;
  $self->{'_regulatory_feature'}  = $regulatory_feature;
  use Scalar::Util qw( weaken );
  # Avoid circular reference and the memory leak that comes with it
  weaken($self->{'_regulatory_feature'});
  return;
}

=head2 regulatory_feature

  Description   : Deprecated, use set_RegulatoryFeature and 
                  get_RegulatoryFeature methods instead.

=cut
sub regulatory_feature {
  my $self               = shift;
  my $regulatory_feature = shift;

  if(defined $regulatory_feature) {
    $self->set_RegulatoryFeature($regulatory_feature);
  }
  return $self->get_RegulatoryFeature;
}

=head2 SO_term

  Example       : print $regulatory_activity->SO_term;

  Description   : Returns the sequence ontology term for this type of 
                  regulatory feature.

  Returntype    : String
  Status        : At risk

=cut

sub SO_term {
  my $self = shift;
  return $self->get_RegulatoryFeature->feature_type->so_accession;
}

sub summary_as_hash {
  my $self = shift;
  
  my $regulatory_feature = $self->get_RegulatoryFeature;
  my $epigenome          = $self->get_Epigenome;
  my $feature_type       = $regulatory_feature->feature_type;

  return {
    regulatory_feature_stable_id => $regulatory_feature->stable_id,
    epigenome           => $epigenome->display_label,
    source              => $regulatory_feature->analysis->logic_name,
    bound_start         => $regulatory_feature->bound_seq_region_start,
    bound_end           => $regulatory_feature->bound_seq_region_end,
    start               => $regulatory_feature->seq_region_start,
    end                 => $regulatory_feature->seq_region_end,
    strand              => $regulatory_feature->strand,
    seq_region_name     => $regulatory_feature->seq_region_name,
    activity            => $self->activity,
    description         => $feature_type->name . ' ' . lc($self->activity) . ' in ' . $epigenome->display_label,
    feature_type        => $feature_type->name,
  };
}

1;
