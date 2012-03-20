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

# Regulatory Features: Cell-specific information
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

__END__

>perl regulatory_features_3.pl

ENSMUSR00000233228:
        4:136390760..136391047-136393044..136393044
        Cell: ESHyb
        Feature Type: Promoter Associated
        Evidence Features:
                Display Label: H3K4me3 - ESHyb Enriched Site; Position: 4:136390760-136391346;
                Display Label: Tcfcp2l1:MA0145.1; Position: 4:136392892-136392905;
ENSMUSR00000233228:
        4:136390453..136391047-136393044..136393044
        Cell: MEF
        Feature Type: Promoter Associated
        Evidence Features:
                Display Label: H3K4me3 - MEF Enriched Site; Position: 4:136390453-136391556;
                Display Label: Tcfcp2l1:MA0145.1; Position: 4:136392892-136392905;
ENSMUSR00000233228:
        4:136391047..136391047-136393044..136393044
        Cell: MultiCell
        Feature Type: Unclassified
        Evidence Features:
                Display Label: DNase1 - ES Enriched Site; Position: 4:136391047-136392800;
                Display Label: Tcfcp2l1 - ES Enriched Site; Position: 4:136392622-136393044;
                Display Label: Tcfcp2l1:MA0145.1; Position: 4:136392892-136392905;
ENSMUSR00000233228:
        4:136391047..136391047-136393044..136393044
        Cell: NPC
        Feature Type: Promoter Associated
        Evidence Features:
                Display Label: H3K4me3 - NPC Enriched Site; Position: 4:136391056-136391320;
                Display Label: Tcfcp2l1:MA0145.1; Position: 4:136392892-136392905;
ENSMUSR00000233228:
        4:136387850..136391047-136393044..136394750
        Cell: ES
        Feature Type: Promoter Associated
        Evidence Features:
                Display Label: H3K27me3 - ES Enriched Site; Position: 4:136387850-136394750;
                Display Label: H3K4me2 - ES Enriched Site; Position: 4:136389904-136390241;
                Display Label: H3K4me3 - ES Enriched Site; Position: 4:136390644-136391578;
                Display Label: H3K4me2 - ES Enriched Site; Position: 4:136390701-136393074;
                Display Label: H3K4me3 - ES Enriched Site; Position: 4:136390797-136392702;
                Display Label: DNase1 - ES Enriched Site; Position: 4:136391047-136392800;
                Display Label: H3K4me3 - ES Enriched Site; Position: 4:136392143-136392603;
                Display Label: Tcfcp2l1 - ES Enriched Site; Position: 4:136392622-136393044;
                Display Label: Tcfcp2l1:MA0145.1; Position: 4:136392892-136392905;
