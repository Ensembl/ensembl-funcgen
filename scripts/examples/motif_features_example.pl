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

Bio::EnsEMBL::Registry->load_registry_from_db(
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
  print $feature->display_label.
     " (".$feature->seq_region_name.":".
    $feature->seq_region_start."-".$feature->seq_region_end.")\n";
}
