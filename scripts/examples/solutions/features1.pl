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
my $gene_adaptor = $efg_db->dnadb->get_GeneAdaptor;
my $regfeat_adaptor = $efg_db->get_RegulatoryFeatureAdaptor;

#Only one gene with this name
my ($gene)  = @{$gene_adaptor->fetch_all_by_external_name('HOXA4')};


my $gene_slice = $gene->feature_Slice->expand('2000', '2000');

my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($gene_slice)};

print 'Found '.scalar(@reg_feats).' RegulatoryFeatures on slice '.$gene_slice->name."\n";

foreach my $reg_feat(@reg_feats){

  print 'Found '.$reg_feat->stable_id.' '.$reg_feat->feature_type->description.' at '.$reg_feat->bound_start.':'.$reg_feat->start.':'.$reg_feat->end.':'.$reg_feat->bound_end."\n";

  foreach my $attr_feat(@{$reg_feat->regulatory_attributes}){
    print 'Attribute Feature '.$attr_feat->feature_type->name."\n";
  }

} 
