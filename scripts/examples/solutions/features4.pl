use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');


#4. We have recently integrated some eQTL data which has enabled the association of Genes and Regulatory features.

#Create a script which retrieves all the Genes for 'chromosome:NCBI36:21:14000000:20000000' and list any which have been linked to RegulatoryFeatures.

#Hint:Use the fetch_all_by_external_name method from the RegulatoryFeatureAdaptor.


my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

my $slice_adaptor = $efg_db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('chromosome', 21, 14000000, 20000000);



#Note how we can use standard adaptor style access for reg feats
#We do have different regulatory FeatureSets, but API defaults to latest.

my $gene_adaptor = $efg_db->dnadb->get_GeneAdaptor;
my $reg_feat_adaptor = $efg_db->get_RegulatoryFeatureAdaptor;


my @genes = @{$gene_adaptor->fetch_all_by_Slice($slice)};

print "got ".scalar(@genes)." genes\n";
my $count = 0;

foreach my $gene(@genes){

  #Grab any RepeatFeatures on this slice
  my @reg_feats = @{$reg_feat_adaptor->fetch_all_by_external_name($gene->stable_id)};

  if(@reg_feats){
	$count ++;
	print "\n".$gene->stable_id.' is in LD with '.scalar(@reg_feats)." RegulatoryFeatures:\n";
  
	#List all the RepeatFeature details
	foreach my $reg_feat(@reg_feats){
	  print $reg_feat->stable_id.':'.$reg_feat->feature_type->name."\n";
	}
  }
}

print "Found $count gene linked Regulatory Features\n";
