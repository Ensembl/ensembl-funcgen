=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::RegulatoryEvidence

=head1 SYNOPSIS
=head1 DESCRIPTION
=head1 SEE ALSO
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

sub evidence_type {
  my $self = shift;
  return $self->_attribute_feature_table;
}

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

1;
