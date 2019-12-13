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

  Bio::EnsEMBL::DBSQL::Funcgen::RegulatoryActivityAdaptor

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

  Bio::EnsEMBL::Registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
      -user => 'anonymous'
  );

  my $regulatory_feature_adaptor       = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryFeature');
  my $regulatory_activity_adaptor      = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryActivity');
  my $regulatory_evidence_link_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryEvidenceLink');

  my $regulatory_feature = $regulatory_feature_adaptor->fetch_by_stable_id('ENSR00000054736');

  print "The regulatory feature with stable id: "  . $regulatory_feature->stable_id . " has the following activities: \n";

  my $regulatory_activity_list = $regulatory_activity_adaptor->fetch_all_by_RegulatoryFeature($regulatory_feature);

  foreach my $current_regulatory_activity (@$regulatory_activity_list) {

    print "The activity in the epigenome "
      . $current_regulatory_activity->get_Epigenome->display_label
      . ' is: '
      . $current_regulatory_activity->activity
      . "\n";

    my $regulatory_evidence_link_list = $regulatory_evidence_link_adaptor->fetch_all_by_RegulatoryActivity($current_regulatory_activity);

    if (@$regulatory_evidence_link_list) {
      print "  It is supported by the following evidence:\n";

      foreach my $current_regulatory_evidence_link (@$regulatory_evidence_link_list) {
        my $evidence = $current_regulatory_evidence_link->get_Evidence;
        print "  - " . $evidence->display_label . ': ' . $evidence->start . '..' . $evidence->end . "\n";
      }
    }
  }

=head1 DESCRIPTION

=cut
package Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryActivityAdaptor;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::RegulatoryActivity';
}

sub _tables {
  return  [ 'regulatory_activity', 'ra' ] ;
}

sub fetch_all_by_regulatory_feature_id {
  my ($self, $regulatory_feature_id) = @_;

  my $constraint = "regulatory_feature_id = $regulatory_feature_id";

  return $self->fetch_all($constraint);
}

sub fetch_all_by_RegulatoryFeature {
  my ($self, $regulatory_feature) = @_;

  return $self->fetch_all_by_regulatory_feature_id($regulatory_feature->dbID);
}

1;
