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

my $feature_set_adaptor = $registry->get_adaptor('Human', 'Funcgen', 'FeatureSet');

# Grab and list all the external feature sets
my @external_feature_sets = @{$feature_set_adaptor->fetch_all_by_feature_class('external')};

foreach my $current_external_feature_set (@external_feature_sets) {
  print "External FeatureSet: " . $current_external_feature_set->name . "\n";
}

my $vista_feature_set = $feature_set_adaptor->fetch_by_name('VISTA enhancer set');

# Now you can get all the features (in this case external features) 
# You can also get features by Slice using get_Features_by_Slice: 
#
foreach my $current_vista_feature (@{$vista_feature_set->get_all_Features}){
    print_feature($current_vista_feature);
    
    # There is no epigenome for these features
    # Feature type indicates vista annotation (eg. active enhancer)
    #
    print $current_vista_feature->feature_type->name."\n";
}

sub print_feature {
  my $feature = shift;
  print $feature->display_label
    . "\t(".$feature->seq_region_name.":"
    . $feature->seq_region_start . "-" . $feature->seq_region_end.")\n";
}
