=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryEvidenceLinkAdaptor

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut
package Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryEvidenceLinkAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Data::Dumper;
use DBI qw(:sql_types);
use Bio::EnsEMBL::Funcgen::RegulatoryEvidenceLink;
use base 'Bio::EnsEMBL::DBSQL::BaseAdaptor';

sub _tables {
  return (
    [ 'regulatory_evidence', 're' ],
  );
}

sub _columns {
  my $self = shift;

  return qw(
    re.regulatory_activity_id
    re.attribute_feature_id
    re.attribute_feature_table
  );
}

sub fetch_all_by_regulatory_activity_id {
  my $self                  = shift;
  my $regulatory_activity_id = shift;

  my $constraint = "re.regulatory_activity_id = ?";
  $self->bind_param_generic_fetch($regulatory_activity_id, SQL_INTEGER);
  
  my $regulatory_evidence_link = $self->generic_fetch($constraint);
  return $regulatory_evidence_link;
}

sub fetch_all_by_RegulatoryActivity {
  my $self               = shift;
  my $regulatory_activity = shift;

  return $self->fetch_all_by_regulatory_activity_id($regulatory_activity->dbID);
}

sub _objs_from_sth {
  my ($self, $sth) = @_;
  
  my (
    $sth_regulatory_activity_id,
    $sth_attribute_feature_id,
    $sth_attribute_feature_table
  );

  $sth->bind_columns (
    \$sth_regulatory_activity_id,
    \$sth_attribute_feature_id,
    \$sth_attribute_feature_table
  );
  
  my @objects_from_sth;
  while ( $sth->fetch() ) {
    my $current_regulatory_evidence_link = Bio::EnsEMBL::Funcgen::RegulatoryEvidenceLink->new;
    
    $current_regulatory_evidence_link->db($self->db);
    $current_regulatory_evidence_link->_regulatory_activity_id($sth_regulatory_activity_id);
    $current_regulatory_evidence_link->_attribute_feature_id($sth_attribute_feature_id);
    $current_regulatory_evidence_link->_attribute_feature_table($sth_attribute_feature_table);
    
    push @objects_from_sth, $current_regulatory_evidence_link;
  }
  return \@objects_from_sth;
}

1;
