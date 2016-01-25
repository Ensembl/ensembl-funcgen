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

my $fset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my $slice_adaptor = $registry->get_adaptor('Human', 'core', 'slice');

my $slice = $slice_adaptor->fetch_by_region('chromosome',1,54960000,54980000);


#Classes of feature set: 'regulatory', 'annotated', 'external'
my @reg_fsets = @{$fset_adaptor->fetch_all_by_type('regulatory')};

foreach my $reg_fset (@reg_fsets) {
	print "Feature Set name: ".$reg_fset->name."\n";
	#Regulatory Feature Sets
	print "\tClass of Feature Set: ".$reg_fset->feature_class."\n";
	#The Regulatory Build
	print "\tAnalysis: ".$reg_fset->analysis->logic_name."\n";
	#Regulatory Feature Type (only makes sense when used together with class)
	print "\tFeature Type: ".$reg_fset->feature_type->name."\n";
	#Regulatory Feature Sets have Cell Type defined
	print "\tCell Type: ".$reg_fset->cell_type->name."\n";
	#Finally, you can also get features from this set
	print "\tSome Features for this Feature Set: \n";
	my @features = @{$reg_fset->get_Features_by_Slice($slice)};
	foreach my $feat (@features) { print "\t\t"; print_feature($feat); }
}



sub print_feature {
	my $feature = shift;
	print 	$feature->display_label. 	
	 	"\t(".$feature->seq_region_name.":".
		$feature->seq_region_start."-".$feature->seq_region_end.")\n";
}
