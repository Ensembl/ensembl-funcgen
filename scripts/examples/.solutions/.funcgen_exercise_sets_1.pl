#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# 1. List the details of all the DataSets in the Human eFG DB
# Create a script which fetches all available DataSets. How many are there? Print out the name, feature type name and cell type name for each.</p>

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
