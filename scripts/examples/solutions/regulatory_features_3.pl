#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

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
