#!/usr/bin/env perl

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

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);


my $slice_adaptor = $registry->get_adaptor('Human', 'Core',    'Slice');
my $slice = $slice_adaptor->fetch_by_region('chromosome', 1, 54_960_000, 54_980_000);

my $regulatory_feature_adaptor = $registry->get_adaptor('Human', 'Funcgen', 'RegulatoryFeature');
my @regulatory_features = @{$regulatory_feature_adaptor->fetch_all_by_Slice($slice)};

foreach my $current_regulatory_feature (@regulatory_features) {
  print $current_regulatory_feature->stable_id.": ";
  print_feature($current_regulatory_feature);
  print "\tFeature Type: ".$current_regulatory_feature->feature_type->name."\n";
}

my $regulatory_feature_demo_stable_id = 'ENSR00000165384';

my $regulatory_activity_adaptor      = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryActivity');
my $regulatory_evidence_link_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryEvidenceLink');

my $regulatory_feature = $regulatory_feature_adaptor->fetch_by_stable_id($regulatory_feature_demo_stable_id);

print "The regulatory feature with stable id: "  . $regulatory_feature->stable_id . " has the following activities: \n";

my $regulatory_activity_list = $regulatory_activity_adaptor->fetch_all_by_RegulatoryFeature($regulatory_feature);

foreach my $current_regulatory_activity (@$regulatory_activity_list) {

  print "\tIn the epigenome "  
    . $current_regulatory_activity->get_Epigenome->display_label 
    . ' it is: ' 
    . $current_regulatory_activity->activity 
    . "\n";

  my $regulatory_evidence_link_list = $regulatory_evidence_link_adaptor->fetch_all_by_RegulatoryActivity($current_regulatory_activity);

  if (@$regulatory_evidence_link_list) {
    print "  This is supported by the following evidence:\n";

    foreach my $current_regulatory_evidence_link (@$regulatory_evidence_link_list) {
      my $evidence = $current_regulatory_evidence_link->get_Evidence;
      print "  - " . $evidence->display_label . ': ' . $evidence->start . '..' . $evidence->end . "\n";
    }
  }
}

sub print_feature {
  my $feature = shift;
  print $feature->display_label.
     " (".$feature->seq_region_name.":".
    $feature->seq_region_start."-".$feature->seq_region_end.")\n";
}
