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
my @motif_features = @{$rf->regulatory_attributes('motif')};

#An example to print motif feature properties
foreach my $motif_feature (@motif_features) {
	print_feature($motif_feature);
	print "\tBinding Matrix: ".$motif_feature->binding_matrix->name."\n";
	print "\tMotif Sequence: ".$motif_feature->seq."\n";
	print "\tMotif Score: ".$motif_feature->score."\n";
	my $afs = $motif_feature->associated_annotated_features();	
	print "\tSupporting Annotated Features\n";
	foreach my $feat (@$afs){
		print "\t\t";
		#Each feature is an annotated feature
		print_feature($feat); 
	}
}

#We can also get motifs via the annotated features.
$rf = $regfeat_adaptor->fetch_by_stable_id('ENSR00000354060'); 
my @annotated_features = @{$rf->regulatory_attributes('annotated')};
foreach my $annot_feat (@annotated_features){
	#One motif feature may appear many times with different annotated features
	my @motif_feats = @{$annot_feat->get_associated_MotifFeatures()};
	if(scalar(@motif_feats)>0){
		print "Annotated Feature with motif feature(s) - ";
		print_feature($annot_feat);
	}
}


sub print_feature {
	my $feature = shift;
	print 	$feature->display_label. 	
	 	"\t(".$feature->seq_region_name.":".
		$feature->seq_region_start."-".$feature->seq_region_end.")\n";
}