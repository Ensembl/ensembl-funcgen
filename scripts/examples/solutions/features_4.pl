#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#1. External Features
#The funcgen DB contains some data directly imported from external sources, like predicted miRNA targets from miRanda, predicted regulatory regions from cisRED, or experimentally verified regulatory regions from EnhancerDB
#Create a script which gets all the 'cisRED motifs' features for the Human region 7:27136000-27139000.
#List their display_labels. List all active regions from EnhancerDB ('VISTA enhancer set') for Human chromosome 7.
#Hint - use their feature type to know which ones are active


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

__END__

>perl features_1.pl

Found 5 cisRED features
craHsap152694
craHsap152679
craHsap152667
craHsap152648
craHsap152626

Found 5 Vista features
hs169
hs170
hs79
hs174
hs183
