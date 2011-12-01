#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $regfeat_adaptor = $registry->get_adaptor('Mouse', 'funcgen', 'regulatoryfeature');

# Regulatory Features: Cell-specific Evidence
# Obtain all the Mouse RegulatoryFeatures associated to the stable ID 'ENSMUSR00000233228' and print their properties
# For each one print the cell type and feature type, and the properties of the evidence features supporting the regulatory feature
# Compare the different cell types. How do you interpret the differences?

my $reg_feature_sets = $regfeat_adaptor->fetch_all_by_stable_ID('ENSMUSR00000233228');
foreach my $rf (@{$reg_feature_sets}){
	print $rf->stable_id.": \n";
	print "\t".$rf->seq_region_name.":".$rf->bound_start."..".$rf->start."-".$rf->end."..".$rf->bound_end."\n";
	print "\tCell: ".$rf->cell_type->name."\n";
	print "\tFeature Type: ".$rf->feature_type->name."\n";
	print "\tEvidence Features: \n";
	map { print_feature($_) } @{$rf->regulatory_attributes()};
}

sub print_feature {
	my $af = shift;
	print "\t\tDisplay Label: ".$af->display_label.";";
	print " Position: ".$af->seq_region_name.":".$af->start."-".$af->end.";";
	print "\n";
}
