use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $fset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my $slice_adaptor = $registry->get_adaptor('Human', 'core', 'slice');

my $slice = $slice_adaptor->fetch_by_region('chromosome',1,54960000,54980000);


#Classes of feature set: 'regulatory', 'annotated', 'external'
my @reg_fsets = @{$fset_adaptor->fetch_all_by_type('regulatory')};

foreach my $reg_fset (@reg_fsets) {
	print "Feature Set name: ".$reg_fset->name."\n";
	#Regulatory Feature Sets
	print "\tClass of Feature Set: ".$reg_fset->feature_class."\n";
	#The Regulatory Build
	print "\tAnalysis: ".$reg_fset->analysis->logic_name."\n";
	#Regulatory Feature Type (only makes sense when used together with class)
	print "\tFeature Type: ".$reg_fset->feature_type->name."\n";
	#Regulatory Feature Sets have Cell Type defined
	print "\tCell Type: ".$reg_fset->cell_type->name."\n";
	#Finally, you can also get features from this set
	print "\tSome Features for this Feature Set: \n";
	my @features = @{$reg_fset->get_Features_by_Slice($slice)};
	foreach my $feat (@features) { print "\t\t"; print_feature($feat); }
}



sub print_feature {
	my $feature = shift;
	print 	$feature->display_label. 	
	 	"\t(".$feature->seq_region_name.":".
		$feature->seq_region_start."-".$feature->seq_region_end.")\n";
}