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

#1. DataSets
#Datasets are a meta container for data, grouping the data obtained by an analysis (FeatureSets) to the underlying raw data (ResultSets).
#Create a script which fetches all available DataSets for Human. How many are there?
#Now get the 'RegulatoryFeatures:MultiCell' data set and print the display label of the product feature set and all the supporting sets.

#Grab the eFG adaptor
my $dsa = $registry->get_adaptor('Human', 'funcgen', 'dataset');

my @human_data_sets = @{$dsa->fetch_all()};
print "There are ".scalar(@human_data_sets)." datasets for Human\n";

#foreach my $dataset (@human_data_sets){
#	print $dataset->name."\t".$dataset->cell_type->name."\t".$dataset->feature_type->name."\n";
#}


my $dset = $dsa->fetch_by_name('RegulatoryFeatures:MultiCell');


#Print the product FeatureSet name
print "Product FeatureSet:\t".$dset->product_FeatureSet->name."\n";

#Now print all the supporting sets
my @ssets = @{$dset->get_supporting_sets};
print "Found ".scalar(@ssets)." supporting sets\n";

foreach my $sset(@ssets){
  print "Supporting ".ucfirst($sset->feature_class)."Feature set:\t".$sset->display_label."\n";

}

__END__

>perl sets_1.pl
There are 612 datasets for Human
Product FeatureSet:     RegulatoryFeatures:MultiCell
Found 245 supporting sets
Supporting AnnotatedFeature set:        DNase1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - NHEK Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - GM06990 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - NHEK Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - K562 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - NHEK Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - IMR90 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HMEC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HSMM Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - NH-A Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - IMR90 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - CD4 Enriched Sites
Supporting AnnotatedFeature set:        Gabp - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cjun - K562 Enriched Sites
Supporting AnnotatedFeature set:        Jund - K562 Enriched Sites
Supporting AnnotatedFeature set:        Max - K562 Enriched Sites
Supporting AnnotatedFeature set:        Nfe2 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Srf - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - K562 Enriched Sites
Supporting AnnotatedFeature set:        Nrsf - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cfos - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - NHEK Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Gabp - GM12878 Enriched Sites
