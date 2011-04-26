use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);


my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');
my $slice_adaptor = $registry->get_adaptor('Human', 'core', 'slice');

my $slice = $slice_adaptor->fetch_by_region('chromosome',1,54960000,54980000);
#Global 'MultiCell' Regulatory Features
my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
foreach my $rf (@reg_feats){ 	
	print $rf->stable_id.": "; 	
	print_feature($rf);
	print "\tCell: ".$rf->cell_type->name."\n"; 	
	print "\tFeature Type: ".$rf->feature_type->name."\n"; 
}

#This gets the global 'MultiCell' Regulatory Feature
print "\nRegulatory Feature ENSR00000165384\n";
my $rf = $regfeat_adaptor->fetch_by_stable_id('ENSR00000165384'); 
#And prints the evidence supporting it
foreach my $feature (@{$rf->regulatory_attributes()}){
	print "\t";
	print_feature($feature);
}

#this gets all cell-specific annotations for this regulatory feature
my $rfs = $regfeat_adaptor->fetch_all_by_stable_ID('ENSR00000165384'); 
foreach my $cell_rf (@{$rfs}){
	#The stable id will always be 'ENSR00000165384' 	
	print $cell_rf->stable_id.": \n"; 	
	#But now it will be for a specific cell type
	print "\tCell: ".$cell_rf->cell_type->name."\n";
	#It will also contain cell-specific annotation
	print "\tType: ".$cell_rf->feature_type->name."\n";
	#And cell-specific extra boundaries
	print 	"\t".$cell_rf->seq_region_name.":".	$cell_rf->bound_start."..".
		$cell_rf->start."-".$cell_rf->end."..".$cell_rf->bound_end."\n";	
	#Unlike the generic MultiCell Regulatory Features, Histone
	# modifications and Polymerase are also used as attributes	
	print "\tEvidence Features: \n"; 	
	foreach my $attr_feat (@{$cell_rf->regulatory_attributes()}){
		print "\t\t";
		print_feature($attr_feat);
	}	
}

# Get all Regulatory Feature Sets
my $fset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my @reg_fsets = @{$fset_adaptor->fetch_all_by_type('regulatory')};
print "Regulatory Feature Sets\n";
foreach my $fset (@reg_fsets) { print "\t".$fset->name."\n"; }

sub print_feature {
	my $feature = shift;
	print 	$feature->display_label. 	
	 	"\t(".$feature->seq_region_name.":".
		$feature->seq_region_start."-".$feature->seq_region_end.")\n";
}

