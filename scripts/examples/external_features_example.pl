#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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
#Grab and list all the external feature sets
my @ext_fsets = @{$fset_adaptor->fetch_all_by_type('external')};

foreach my $ext_fset (@ext_fsets){
  print "External FeatureSet: ".$ext_fset->name."\n";
}

#Grab the specific Vista set
my $vista_fset = $fset_adaptor->fetch_by_name('VISTA enhancer set');

#Now you can get all the features (in this case external features) 
#You can also get features by Slice using get_Features_by_Slice: 
foreach my $vista_feature (@{$vista_fset->get_all_Features()}){
	print_feature($vista_feature);
	#There is no cell type for these features
	#Feature type indicates vista annotation (eg. active enhancer)
	print "\tFeature Type: ".$vista_feature->feature_type->name."\n";
}

sub print_feature {
	my $feature = shift;
	print 	$feature->display_label. 	
	 	"\t(".$feature->seq_region_name.":".
		$feature->seq_region_start."-".$feature->seq_region_end.")\n";
}
