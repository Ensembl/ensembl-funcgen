#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $rfa = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');

# Motif features represent putative binding sites. They are always associated to annotated features.
# Get the 'motif' regulatory attributes associated to the human Regulatory Feature 'ENSR00001227187'. Print their properties 

my $rf = $rfa->fetch_by_stable_id('ENSR00001227187');
foreach my $mf ( @{$rf->regulatory_attributes('motif')}){
	print "Display Label: ".$mf->display_label.";";
	print " Position: ".$mf->seq_region_name.":".$mf->start."-".$mf->end.";";
	print " Score: ".$mf->score.";";
	print "\n";
	
	# Print the properties of the annotated features associated to each motif feature.
	# Print the name of the feature set associated to these annotated features.
	foreach my $af (@{$mf->associated_annotated_features()}){
		print "\tAssociated Feature Set: ".$af->feature_set->name."\n";
		print_feature($af);		
	}
	
};

sub print_feature {
	my $af = shift;
	print "\tCell: ".(defined($af->cell_type) ? $af->cell_type->name : "Undetermined").";"; 
	print " Feature Type: ".$af->feature_type->name.";";
	print " Position: ".$af->seq_region_name.":".$af->start."-".$af->end.";";
	print " Score: ".$af->score.";";
	print " Summit: ".$af->summit.";";
	print "\n";	
}
