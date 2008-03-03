use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

#Investigate what DataSets are present in the data base.
#Create a script which lists all the DataSet names, feature types, cell types and their supporting sets


my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

my $slice_adaptor = $efg_db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_name('chromosome:NCBI36:7:27132000:27139000');

my $featureset_adaptor = $efg_db->get_FeatureSetAdaptor;
my $cisred_fset = $featureset_adaptor->fetch_by_name('cisRED group motifs');


my @feats = @{$cisred_fset->get_Features_by_Slice($slice)};

print 'Found '.scalar(@feats).' '.$cisred_fset->name.' on slice '.$slice->name."\n";



foreach my $feat(@feats){

  print 'Found '.$feat->display_label."\n";

} 
