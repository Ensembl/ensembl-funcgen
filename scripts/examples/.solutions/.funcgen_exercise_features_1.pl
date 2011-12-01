#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#2. External Features
# The eFG contains some data directly imported from external sources, like predicted miRNA targets from miRanda, predicted regulatory regions from cisRED, experimentally verified regulatory regions from enhancerDB
#Create a script which gets all the 'cisRED motifs' features for the region 7:27136000-27139000. List their display_labels.
# Each External Feature can have their own specific feature type (like an extra attribute)
#List all active regions from enhancerDB ('VISTA enhancer set') for chromosome 7 (use the feature type to see which ones have an Enhancer activity)

#Grab the eFG adaptor
my $sa = $registry->get_adaptor('Human', 'core', 'slice');
my $fsa = $registry->get_adaptor('Human', 'funcgen', 'featureset');

my $cisred_slice = $sa->fetch_by_region('chromosome','7',27136000,27139000);
my $cisred_fset = $fsa->fetch_by_name('cisRED motifs');
my @cisred_afs = @{$cisred_fset->get_Features_by_Slice($cisred_slice)};
print "Found ".scalar(@cisred_afs)." cisRED features\n";
foreach my $af (@cisred_afs){
	print $af->display_label."\n";
}


my $vista_fset = $fsa->fetch_by_name('VISTA enhancer set');
my $vista_slice = $sa->fetch_by_region('chromosome','7');
my @vista_afs = @{$vista_fset->get_Features_by_Slice($vista_slice)};

print "Found ".scalar(@vista_afs)." Vista features\n";
foreach my $af (@vista_afs){
	if($af->feature_type->name =~ /Enhancer/){ print $af->display_label."\n"; }
}
