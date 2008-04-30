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
