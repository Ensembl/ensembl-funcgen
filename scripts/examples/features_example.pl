use strict;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db
  (
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
  );


#Grab the adaptors and a slice to play with
my $efg_db             = $registry->get_DBAdaptor('Human', 'funcgen');
my $featureset_adaptor = $efg_db->get_FeatureSetAdaptor;
my $slice_adaptor      = $efg_db->get_SliceAdaptor;
my $slice              = $slice_adaptor->fetch_by_region('chromosome', '17', 570000, 600000);

#list the features for the external fset 'miRanda miRNA'
my $fset     = $featureset_adaptor->fetch_by_name('miRanda miRNA');
my @features = @{$fset->get_Features_by_Slice($slice)};

foreach my $feat(@features){
  print $feat->display_label."\t".$feat->feature_Slice->name."\n";
}


#Get some RegulatoryFeatures and list their attributes and feature set name
my $regfeat_adaptor = $efg_db->get_RegulatoryFeatureAdaptor;
my @reg_features    = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};


#Print out some info
foreach my $reg_feat(@reg_features){
  
  print "\n".$reg_feat->stable_id.' '.$reg_feat->feature_type->name."\n";

  foreach my $attr_feat(@{$reg_feat->regulatory_attributes}){
	print 'AttributeFeature '.$attr_feat->feature_type->name."\n";
  }
}
