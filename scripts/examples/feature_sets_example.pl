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

my $feature_set_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my $slice_adaptor = $registry->get_adaptor('Human', 'Core', 'Slice');
my $slice = $slice_adaptor->fetch_by_region('chromosome', 1, 54_960_000, 54_980_000);

my $annotated_feature_sets = $feature_set_adaptor->fetch_all_by_feature_class('annotated');

foreach my $current_annotated_feature_set (@$annotated_feature_sets) {
    print "---- " . $current_annotated_feature_set->name . " ----\n";
    print $current_annotated_feature_set->feature_class."\n";
    print $current_annotated_feature_set->analysis->logic_name."\n";
    print $current_annotated_feature_set->feature_type->name."\n";
    print $current_annotated_feature_set->epigenome->name."\n";
    
    # Finally, you can also get features from this set
    #
    my @annotated_features = @{$current_annotated_feature_set->get_Features_by_Slice($slice)};
    foreach my $current_feature (@annotated_features) {
      print_feature($current_feature);
    }
    print "\n";
}

sub print_feature {
  my $feature = shift;
  print "\t" . $feature->display_label
    . "\t(".$feature->seq_region_name.":"
    . $feature->seq_region_start . "-" . $feature->seq_region_end.")\n";
}
