#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# 2. FeatureSets
# Feature Sets hold processed data.
# Print the name of the feature sets for the 'GM12878' cell type.
# Print the name of the feature sets for the 'CTCF' feature type.
# Is the FeatureSet 'VISTA enhancer set' associated to any cell type of feature type? 
# Is there a DataSet for the FeatureSet 'VISTA enhancer set'? Hint: Use fetch_by_product_FeatureSet. Any idea why this is so?
# Is there a DataSet for the FeatureSet 'RegulatoryFeatures:MultiCell'? Any idea why this is so?

#Grab the eFG adaptor
my $fsa = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my $cta = $registry->get_adaptor('Human', 'funcgen', 'celltype');


#Print GM12878 cell type details
my @GM12878_feature_sets = @{$fsa->fetch_all_by_CellType($cta->fetch_by_name('GM12878'))};
print "There are ".scalar(@GM12878_feature_sets)." feature sets for GM12878:\n";
foreach my $featureset (@GM12878_feature_sets){
	print "\t".$featureset->name."\n";
}

#Print CTCF feature type details
my $fta = $registry->get_adaptor('Human', 'funcgen', 'featuretype');
my @CTCF_feature_sets = @{$fsa->fetch_all_by_FeatureType($fta->fetch_by_name('CTCF'))};
print "\n\nThere are ".scalar(@CTCF_feature_sets)." feature sets for CTCF:\n";
foreach my $featureset (@CTCF_feature_sets){
	print "\t".$featureset->name."\n";
}


#Print VISTA sets details
my $vista_set = $fsa->fetch_by_name('VISTA enhancer set');
print "\n\n".$vista_set->display_label."\n";

my $cell_type = ($vista_set->cell_type) ? $vista_set->cell_type->name : 'NON-DEFINED';
print "\tCell type:\t".$cell_type;

my $ftype = ($vista_set->feature_type) ? $vista_set->feature_type->name : 'NON-DEFINED';
print "\tFeatureType:\t$ftype\n";

my $dsa     = $registry->get_adaptor('Human', 'funcgen', 'dataset');
my $dset    = $dsa->fetch_by_product_FeatureSet($vista_set);
my $dataset = ($dset) ? $dset->name : 'NON-DEFINED';
print "\tDataSet:\t$dataset\n";

#Print RegFeat details
my $multicellreg_set = $fsa->fetch_by_name('RegulatoryFeatures:MultiCell');
print "\n\n".$multicellreg_set->display_label."\n";
$dset    = $dsa->fetch_by_product_FeatureSet($multicellreg_set);
$dataset = ($dset) ? $dset->name : 'NON-DEFINED';
print "\tDataSet:\t".$dataset."\n";
