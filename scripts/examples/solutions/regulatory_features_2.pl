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



# Regulatory Features: What RegulatoryFeatures are near the oncogene BRCA2?
# Create a script which gets all the RegulatoryFeatures within 1KB of the 'BRCA2' gene.
# Print out their stable IDs, bound_start/end and start/end values, name of the cell and feature types.


# HINT: Use fetch_all_by_external_name with 'BRCA2' to get the gene object.
# HINT: Look at the argument for fetch_by_gene_stable_id or use the Gene->feature_Slice and Slice->expand methods.


#Get the Gene and RegualtoryFeature adaptors
my $gene_adaptor = $registry->get_adaptor('Human', 'core', 'gene');
my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');

my $gene_name = 'BRCA2';

my @genes = @{$gene_adaptor->fetch_all_by_external_name($gene_name)};
print scalar(@genes)." human gene(s) named $gene_name\n";

my $gene = $genes[1];
my $slice = $registry->get_adaptor('human', 'core', 'slice')->fetch_by_gene_stable_id($gene->stable_id, 1000);

#my $slice = $b->feature_Slice;
#print "Slice corresponding to $gene_name: ".$slice->seq_region_name."\t".$slice->start."\t".$slice->end."\t".$slice->strand."\n";
#$slice = $slice->expand(1000,1000);

print "Slice 1Kb around $gene_name: ".$slice->seq_region_name."\t".$slice->start."\t".$slice->end."\t".$slice->strand."\n";

my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
print scalar(@reg_feats)." Regulatory Features 1kb around $gene_name\n";


foreach my $rf (@reg_feats){
	print $rf->stable_id.": \n";
	print "\t".$rf->seq_region_name.":".$rf->bound_start."..".$rf->start."-".$rf->end."..".$rf->bound_end."\n";
	print "\tCell: ".$rf->cell_type->name."\n";
	print "\tFeature Type: ".$rf->feature_type->name."\n";
	#print "\tEvidence Features: \n";

	#map { print_feature($_) } @{$rf->regulatory_attributes()};
}

#sub print_feature {
#	my $af = shift;
#	print "\t\tDisplay Label: ".$af->display_label.";";
#	print " Position: ".$af->seq_region_name.":".$af->start."-".$af->end.";";
#	print "\n";
#}

__END__

>perl regulatory_features_2.pl
1 human gene(s) named BRCA2
Slice 1Kb around BRCA2: 13      32888611        32974805        1
11 Regulatory Features 1kb around BRCA2
ENSR00000054736: 
        13:370..370-1921..1921
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001036089: 
        13:3636..3636-3911..3911
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001508453: 
        13:6120..6120-6449..6449
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00000513985: 
        13:9391..9391-9889..9889
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001508454: 
        13:15917..15917-16386..16386
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001036090: 
        13:39644..39644-39987..39987
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00000274472: 
        13:50507..50507-50664..50664
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001508455: 
        13:53399..53399-53568..53568
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001508456: 
        13:67051..67051-67514..67514
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001036091: 
        13:70251..70251-70573..70573
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00000513988: 
        13:76809..76809-77091..77091
        Cell: MultiCell
        Feature Type: Unclassified
