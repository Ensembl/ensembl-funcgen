#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
