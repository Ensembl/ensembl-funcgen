use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

#What are the supporting sets for the 'RegulatoryFeatures' DataSet?
#Create a script which fetches the 'RegulatoryFeatures' DataSet. Print out the name, feature type and cell type of each of the supporting sets.


my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

my $dataset_adaptor = $efg_db->get_DataSetAdaptor;

my $data_set = $dataset_adaptor->fetch_by_name('RegulatoryFeatures');

print $dset->name.' contains the following '.$dset->supporting_set_type." feature supporting sets:\n";

foreach my $sset(@{$data_set->get_supporting_sets}){
  print "Supporting set:\t".$sset->name."\t".ucfirst($sset->type)."Features\n";
} 

my $product_fset = $data_set->product_FeatureSet;

print 'The product FeatureSet for this DataSet is '.$product_fset->name.' containing '.ucfirst($product_fset->type)."Features\n";
