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

my $rfa = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');

#3. Motif Features
#Motif features represent putative binding sites based on alignments of PWMs from JASPAR. MotifFeatures are always associated to AnnotatedFeatures representing Transcription Factor (TF) Binding. More information about how we integrate these into the regulatory build process can be found here.
#Get the 'motif' regulatory attributes associated to the Human Regulatory Feature 'ENSR00001227187'. Print their properties.
#Hint: use 'motif' as a parameter for regulatory_attributes.
#Print the properties of the annotated features associated to the motif feature.

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

__END__


>perl features_3.pl

Display Label: CTCF:MA0139.1; Position: 6:132128805-132128823; Score: 0.863;
        Associated Feature Set: HUVEC_CTCF_ENCODE_Uw_SWEMBL_R015
        Cell: HUVEC; Feature Type: CTCF; Position: 6:132128612-132128993; Score: 22.126581; Summit: 132128815;
        Associated Feature Set: H1ESC_CTCF_ENCODE_Uta_SWEMBL_R015
        Cell: H1ESC; Feature Type: CTCF; Position: 6:132128668-132129004; Score: 50.586062; Summit: 132128831;
        Associated Feature Set: HepG2_CTCF_ENCODE_Uta_SWEmbl_R015_D150
        Cell: HepG2; Feature Type: CTCF; Position: 6:132128676-132128987; Score: 27.973154; Summit: 132128824;
        Associated Feature Set: HepG2_CTCF_ENCODE_Uw_SWEmbl_R015_D150
        Cell: HepG2; Feature Type: CTCF; Position: 6:132128678-132129034; Score: 82.49168; Summit: 132128828;

Display Label: Nrsf:MA0138.1; Position: 6:132129499-132129517; Score: 0.805;
        Associated Feature Set: H1ESC_Nrsf_ENCODE_Hudsonalpha_SWEMBL_R015
        Cell: H1ESC; Feature Type: Nrsf; Position: 6:132129342-132129713; Score: 103.256328; Summit: 132129519;
