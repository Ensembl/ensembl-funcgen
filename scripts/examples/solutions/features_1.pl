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

#Annotated Features
#Annotated Features represents the results of an analysis of raw or processing signal data. These correspond to regions in the genome enriched for specific events (like TF binding or Histone Marks) i.e. they are 'peak calls'.
#
#Fetch the AnnotatedFeatures in the region Y:5000000-40000000 for the Human feature sets:
#    K562_DNase1_ENCODE_Duke_SWEmbl_R0025_D150<br>
#    HepG2_DNase1_ENCODE_Duke_SWEmbl_R0025_D150
#
#Print the number of features returned by each, the details of the CellType associated with the FeatureSet and the details of the few features including the 'summit'.
#
#What are the differences and why?
#
#>Optional: Print the properties(logic_name, display_label, parameters) of the <a href="http://www.ensembl.org/info/docs/Doxygen/core-api//classBio_1_1EnsEMBL_1_1Analysis.html">Analysis</a> used in the previous feature sets.



#Grab the eFG adaptor
my $sa = $registry->get_adaptor('Human', 'core', 'slice');
my $fsa = $registry->get_adaptor('Human', 'funcgen', 'featureset');

# You can get it by the name
#Carefull: some parts of the Y chomosome are PARs i.e. same as the X...
my $slice = $sa->fetch_by_region('chromosome','Y',5000000,40000000);

# Compare the number of annotated features in the region Y:5000000-40000000 between the Human feature sets 'K562_DNase1_ENCODE_Duke_SWEmbl_R0025_D150' and 'HepG2_DNase1_ENCODE_Duke_SWEmbl_R0025_D150'.Investigate the differences (hint: check cell types).
my $K562_DNAse_fset = $fsa->fetch_by_name('K562_DNase1_ENCODE_Duke_SWEmbl_R0025_D150');
my @K562_DNAse_afs = @{$K562_DNAse_fset->get_Features_by_Slice($slice)};

print "K562 contains ".scalar(@K562_DNAse_afs)." features\n";
print $K562_DNAse_fset->cell_type->name." has gender ".$K562_DNAse_fset->cell_type->gender."\n";



my $HepG2_DNAse_fset = $fsa->fetch_by_name('HepG2_DNase1_ENCODE_Duke_SWEmbl_R0025_D150');
my @HepG2_DNAse_afs = @{$HepG2_DNAse_fset->get_Features_by_Slice($slice)};

print "HepG2 contains ".scalar(@HepG2_DNAse_afs)." features\n";
print $HepG2_DNAse_fset->cell_type->name." has gender ".$HepG2_DNAse_fset->cell_type->gender."\n";

my $count = 0;

foreach my $af(@HepG2_DNAse_afs){
	print "Start: ".$af->seq_region_start." End: ".$af->seq_region_end." Score: ".$af->score." Summit: ".$af->summit."\n";
  $count ++;
  last if $count == 5;
}

print "\n";



# (optional) Print the properties of the analysis used in the previous feature sets
print "K562 DNAse\n";
print_analysis($K562_DNAse_fset->analysis);
print "\n";

print "HepG2 DNAse\n";
print_analysis($HepG2_DNAse_fset->analysis);
print "\n";


sub print_analysis {
	my $an = shift;
	print "Logic Name: ".$an->logic_name."\n";
	print "Display Label: ".$an->display_label."\n";
	print "Program: ".$an->program."\t".$an->parameters."\n";
}

__END__


>perl features_1.pl
K562 contains 0 features
K562 has gender female
HepG2 contains 54 features
HepG2 has gender male
Start: 5505833 End: 5506015 Score: 33.063714 Summit: 5505924
Start: 6778390 End: 6778689 Score: 28.400612 Summit: 6778540
Start: 7141665 End: 7142190 Score: 49.68536 Summit: 7141812
Start: 9929940 End: 9930295 Score: 23.551168 Summit: 9930035
Start: 9930760 End: 9930975 Score: 37.184378 Summit: 9930867

K562 DNAse
Logic Name: SWEmbl_R0025_D150
Display Label: SWEmbl
Program: SWEmbl -f 150 -R 0.0025 -d 150

HepG2 DNAse
Logic Name: SWEmbl_R0025_D150
Display Label: SWEmbl
Program: SWEmbl -f 150 -R 0.0025 -d 150
