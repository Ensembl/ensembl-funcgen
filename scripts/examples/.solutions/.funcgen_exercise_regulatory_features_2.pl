#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');

# Regulatory Features: Evidence
# Print out the display_label, start/end values of all the evidence features for Human regulatory feature "ENSR00000623613".
# Compare with the start/end values of the regulatory feature itself.
my $rf = $regfeat_adaptor->fetch_by_stable_id('ENSR00000623613');
print $rf->stable_id.": \n";
print "\t".$rf->seq_region_name.":".$rf->bound_start."..".$rf->start."-".$rf->end."..".$rf->bound_end."\n";
print "\tEvidence Features: \n";
map { print_feature($_) } @{$rf->regulatory_attributes()};

sub print_feature {
	my $af = shift;
	print "\t\tDisplay Label: ".$af->display_label.";";
	print " Position: ".$af->seq_region_name.":".$af->start."-".$af->end.";";
	print "\n";
}
