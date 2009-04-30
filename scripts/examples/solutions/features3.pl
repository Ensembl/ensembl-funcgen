use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

#3. Generate a Promoter profile
# Fetch all the promoter associated regulatory features for chr??
# Count the presence of each supporting feature type to identify what defines a promoter

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

my $slice_adaptor = $efg_db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_name('chromosome:NCBI36:1:500000:1000000');

my $featureset_adaptor = $efg_db->get_FeatureSetAdaptor;
my $regfeat_fset = $featureset_adaptor->fetch_by_name('RegulatoryFeatures');


my @reg_feats = @{$regfeat_fset->get_Features_by_Slice($slice)};

my %feature_type_counts;


foreach my $reg_feat(@reg_feats){
  my %seen_feature_types;

  next if $reg_feat->feature_type->name !~ /Promoter/;
  
  #Get all supporting features for this reg feat
  foreach my $reg_attr_feature(@{$reg_feat->regulatory_attributes}){

	#Only count if we have not already seen it.
	if(! exists $seen_feature_types{$reg_attr_feature->feature_type->name}){
	  $feature_type_counts{$reg_attr_feature->feature_type->name}++;
	}

  }

} 

print 'Promoter profile from '.scalar(@reg_feats)." RegulatoryFeatures:\n";


foreach my $ftype(keys(%feature_type_counts)){
  printf("$ftype% 11d\n", $feature_type_counts{$ftype});
}
