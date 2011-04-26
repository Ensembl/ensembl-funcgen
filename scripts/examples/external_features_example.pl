use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $fset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');
#Grab and list all the external feature sets
my @ext_fsets = @{$fset_adaptor->fetch_all_by_type('external')};

foreach my $ext_fset (@ext_fsets){
  print "External FeatureSet: ".$ext_fset->name."\n";
}

#Grab the specific Vista set
my $vista_fset = $fset_adaptor->fetch_by_name('VISTA enhancer set');

#Now you can get all the features (in this case external features) 
#You can also get features by Slice using get_Features_by_Slice: 
foreach my $vista_feature (@{$vista_fset->get_all_Features()}){
	print_feature($vista_feature);
	#There is no cell type for these features
	#Feature type indicates vista annotation (eg. active enhancer)
	print "\tFeature Type: ".$vista_feature->feature_type->name."\n";
}

sub print_feature {
	my $feature = shift;
	print 	$feature->display_label. 	
	 	"\t(".$feature->seq_region_name.":".
		$feature->seq_region_start."-".$feature->seq_region_end.")\n";
}
