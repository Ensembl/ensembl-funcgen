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
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db
  (
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
  );


#Grab the adaptors
my $efg_db           = $registry->get_DBAdaptor('Human', 'funcgen');
my $dataset_adaptor  = $efg_db->get_DataSetAdaptor;

#Grab the CTCF Nessie data set
my $data_set         = $dataset_adaptor->fetch_by_name('Nessie_NG_STD_2_ctcf_ren_BR1');

#Grab underlying supporting sets
my @supporting_sets  = @{$data_set->get_supporting_sets};

#Print some info about the supporting sets and feature set
foreach my $sset ( @supporting_sets ){
  print "Supporting set:\t".$sset->name."\n";
  print 'Produced by analysis '.$sset->analysis->logic_name."\n";
}

my $pfset = $data_set->product_FeatureSet;

print 'Product FeatureSet is '.$pfset->name."\n";
print 'Produced by analysis '.$pfset->analysis->logic_name."\n";


#Grab and list all the external feature sets
my $featureset_adaptor = $efg_db->get_FeatureSetAdaptor;
my @ext_fsets = @{$featureset_adaptor->fetch_all_by_type('external')};

foreach my $ext_fset(@ext_fsets){
  print "External FeatureSet:\t".$ext_fset->name."\n";
}
