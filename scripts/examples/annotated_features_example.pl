use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');

my $rf = $regfeat_adaptor->fetch_by_stable_id('ENSR00000165384'); 
my @annotated_features = @{$rf->regulatory_attributes('annotated')};

#An example to print annotated feature properties
foreach my $annotated_feature (@annotated_features) {
	print_feature($annotated_feature);
	print "\tFeature Type: ".$annotated_feature->feature_type->name."\n";
	print "\tFeature Set: ".$annotated_feature->feature_set->name."\n";
	#Analysis-depends property
	print "\tScore: ".$annotated_feature->score."\n";
	#Summit is usually present, but may not exist
	print "\tSummit: ".$annotated_feature->summit."\n" if defined($annotated_feature->summit);
}

#Get all Annotated Feature Sets
my $fset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my @annot_fsets = @{$fset_adaptor->fetch_all_by_type('annotated')};
print "There are ".scalar(@annot_fsets)." Annotated Feature Sets\n";

sub print_feature {
	my $feature = shift;
	print 	$feature->display_label. 	
	 	"\t(".$feature->seq_region_name.":".
		$feature->seq_region_start."-".$feature->seq_region_end.")\n";
}