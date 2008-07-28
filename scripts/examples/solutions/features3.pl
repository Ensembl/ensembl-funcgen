use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');


#3. Recent research has shown that repeat regions can induce secondary structure which alters promoter viability?

#Create a script which retrieves all RegulatoryFeatures for the slice 'chromosome:NCBI36:17:5000000:7000000'. Identify promoter associated regions which contain repeat regions and list the repeat types.

#Hint: Promoter FeatureType names are 'Promoter Associated' and 'Promoter Associated - Cell type specific'



my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

my $slice_adaptor = $efg_db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('chromosome', '17', 5000000, 7000000);


#Note how we can use standard adaptor style access for reg feats
#We do have different regulatory FeatureSets, but API defaults to latest.

my $rep_feat_adaptor = $efg_db->dnadb->get_RepeatFeatureAdaptor;
my $reg_feat_adaptor = $efg_db->get_RegulatoryFeatureAdaptor;
my @reg_feats = @{$reg_feat_adaptor->fetch_all_by_Slice($slice)};

print 'Found '.scalar(@reg_feats).' on slice '.$slice->name."\n";

my $rep_count = 0;

foreach my $reg_feat(@reg_feats){

  next if $reg_feat->feature_type->name !~ /Promoter Associated/;

  #Grab the slice of the RegulatoryFeature core region.
  my $feat_slice = $reg_feat->feature_Slice;

  #Grab any RepeatFeatures on this slice
  my @rep_feats = @{$rep_feat_adaptor->fetch_all_by_Slice($feat_slice)};

  if(@rep_feats){
	$rep_count ++;
	print "\n".$reg_feat->stable_id.' contains '.scalar(@rep_feats)." RepeatFeatures:\n";
  
	#List all the RepeatFeature details
	foreach my $rep_feat(@rep_feats){

	  my $rep_cons = $rep_feat->repeat_consensus;
	  print $rep_cons->name.':'.$rep_cons->repeat_class.':'.$rep_cons->repeat_type.':'.$rep_cons->repeat_consensus."\n";
	}
  }
}

print "Found $rep_count Promoter Associated Regulatory Features containing repeats\n";
