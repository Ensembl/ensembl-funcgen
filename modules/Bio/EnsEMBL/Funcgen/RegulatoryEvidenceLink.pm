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

  Bio::EnsEMBL::Funcgen::RegulatoryEvidenceLink

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
  );

  my $regulatory_evidence_link_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryEvidenceLink');
  my $regulatory_evidence_link_list = $regulatory_evidence_link_adaptor->fetch_all_by_regulatory_activity_id(3);

  foreach my $current_regulatory_evidence_link (@$regulatory_evidence_link_list) {
    my $evidence = $current_regulatory_evidence_link->get_Evidence;
    print $evidence->display_label . ': ' . $evidence->start . '..' . $evidence->end . "\n";

  }

=head1 DESCRIPTION

  Object to handle evidence for the activities of regulatory features.

=cut

package Bio::EnsEMBL::Funcgen::RegulatoryEvidenceLink;

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

sub _regulatory_activity_id {
  my $self = shift;
  my $regulatory_activity_id = shift;

  if ($regulatory_activity_id) {
    $self->{'_regulatory_activity_id'}  = $regulatory_activity_id;
  }
  return $self->{'_regulatory_activity_id'};
}

sub _attribute_feature_id {
  my $self = shift;
  my $attribute_feature_id = shift;

  if ($attribute_feature_id) {
    $self->{'_attribute_feature_id'}  = $attribute_feature_id;
  }
  return $self->{'_attribute_feature_id'};
}

sub _attribute_feature_table {
  my $self = shift;
  my $attribute_feature_table = shift;

  if ($attribute_feature_table) {
    $self->{'_attribute_feature_table'}  = $attribute_feature_table;
  }
  return $self->{'_attribute_feature_table'};
}

=head2 evidence_type

  Example    : if ($regulatory_evidence_link->evidence_type eq 'motif') { 
                 print "Doing something with motif evidence.\n"; 
               };

  Description: Returns either 'motif' or 'annotated', depending on whether 
               the evidence is linking to Bio::EnsEMBL::Funcgen::MotifFeature 
               or Bio::EnsEMBL::Funcgen::AnnotatedFeature

  Returntype : String, either 'motif' or 'annotated'.
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub evidence_type {
  my $self = shift;
  return $self->_attribute_feature_table;
}

=head2 get_Evidence

  Example    : my $evidence = $regulatory_evidence_link->get_Evidence;

  Description: Returns the object this is linking to.

  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature or 
               Bio::EnsEMBL::Funcgen::AnnotatedFeature
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub get_Evidence {
  my $self = shift;

  if($self->evidence_type eq 'motif') {
    return $self
      ->db
      ->get_MotifFeatureAdaptor
      ->fetch_by_dbID(
	$self->_attribute_feature_id
      );
  }
  if($self->evidence_type eq 'annotated') {
    return $self
      ->db
      ->get_AnnotatedFeatureAdaptor
      ->fetch_by_dbID(
	$self->_attribute_feature_id
      );
  }
  throw("Unknown evidence type: " . $self->evidence_type);
}

=head2 get_Evidence_on_Slice

  Arg 1      : (mandatory) Bio::EnsEMBL::Slice

  Example    : my $regulatory_feature_slice = $regulatory_feature->slice;
               my $regulatory_evidence_link = $self->get_RegulatoryEvidenceLink;
               my @regulatory_evidence = map {
                 $_->get_Evidence_on_Slice($regulatory_feature_slice) 
               } @$regulatory_evidence_link;

  Description: Returns the object this is linking to on a given slice.
  
               The difference to get_Evidence is that this method accepts a 
               slice object and the evidence returned will be on the same 
               slice.

               This is useful, if you want to make sure that the coordinates 
               match up. E.g.: When you want to use the coordinates of the
               motif feature evidence to find bases in a regulatory feature.

  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature or 
               Bio::EnsEMBL::Funcgen::AnnotatedFeature
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub get_Evidence_on_Slice {
  my $self  = shift;
  my $slice = shift;

  if($self->evidence_type eq 'motif') {
    return $self
      ->db
      ->get_MotifFeatureAdaptor
      ->fetch_all_by_Slice_constraint(
	$slice,
	'motif_feature_id=' . $self->_attribute_feature_id
	# The constraint is an id, so this will return an array with one
	# element. Dereferencing here, so this returns the object only.
      )->[0];
  }
  if($self->evidence_type eq 'annotated') {
    return $self
      ->db
      ->get_AnnotatedFeatureAdaptor
      ->fetch_all_by_Slice_constraint(
	$slice,
	'annotated_feature_id=' . $self->_attribute_feature_id
      )->[0];
  }
  throw("Unknown evidence type: " . $self->evidence_type);
}

1;
