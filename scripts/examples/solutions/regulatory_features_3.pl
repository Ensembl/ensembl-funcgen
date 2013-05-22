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

# Regulatory Features: associated expermental evidence
# Now for just the ENSR00000623613 MultiCell feature, print out the display_label, start/end values of all the evidence features
# Compare with the start/end values of the regulatory feature itself.
# HINT: By default the fetch_by_stable_id method returns just the MultiCell features.

my $rf = $regfeat_adaptor->fetch_by_stable_id('ENSR00000623613');

print $rf->stable_id.":\t\t\t\t\t\t".
  $rf->seq_region_name.":".$rf->bound_start."..".$rf->start."-".$rf->end."..".$rf->bound_end."\n";
print "\tEvidence Features: \n";

my @attributes = @{$rf->regulatory_attributes()};

print "Found ".scalar(@attributes)."\n";

map { print_feature($_) } @attributes;

sub print_feature {
  my $af = shift;
  print "\t\tDisplay Label: ".$af->display_label.";";
  print " Position:\t".$af->seq_region_name.":".$af->start."-".$af->end.";";
  print "\n";
}


__END__

>perl regulatory_features_3.pl
