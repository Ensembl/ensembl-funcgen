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

#4. Binding Matrices and motif strength
#Each MotifFeature is associated with a PWM, which are represented by the 'BindingMatrix' class. The MotifFeature score represents the relative binding affinity with respect to the PWM defined in the BindingMatrix.
#Using the Motif feature obtained in exercise 3, get the associated Binding Matrix and print some details.
#Check potential effect of changes in the sequence of the motif feature on the relative strength of that motif feature.
#Check the GERP conservation scores along the motif. Compare with the JASPAR matrix.

my $rf = $rfa->fetch_by_stable_id('ENSR00001227187');
my $mf = $rf->regulatory_attributes('motif')->[0];
my $slice = $mf->feature_Slice;
my $seq = $mf->seq;

my $bm = $mf->binding_matrix;
print "Relative affinity of the ".$bm->name." motif\t$seq :".$bm->relative_affinity($seq)."\n";
print "Relative affinity of changed motif\tTGGTGTGCAGAGGGCATGA :".$bm->relative_affinity("TGGTGTGCAGAGGGCATGA")."\n";
print "                                            ^\n";
print "Relative affinity of changed motif\tTGGTCTGCGGAGGGCATGA :".$bm->relative_affinity("TGGTCTGCGGAGGGCATGA")."\n";
print "                                                ^\n";




#http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?ID=MA0139.1&rm=present&collection=CORE


#get method_link_species_set object for gerp conservation scores for mammals
my $mlss_adaptor = $registry->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");

#Get cons scores for given slice
my $display_size = $slice->end - $slice->start + 1; 
my $cs_adaptor = $registry->get_adaptor("Multi", 'compara', 'ConservationScore');
my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $display_size);
print join("\t",split('',$seq))."\n";


map { printf("%.3f\t", $_->diff_score); } @$scores;
print "\n";
print "                                ^\n";

print $bm->frequencies;
print "                  ^\n";
print "                                  ^\n";

__END__


>perl features_4.pl

Relative affinity of the MA0139.1 motif TGGTCTGCAGAGGGCATGA :0.863462959414993
Relative affinity of changed motif      TGGTGTGCAGAGGGCATGA :0.747783950261183
                                            ^
Relative affinity of changed motif      TGGTCTGCGGAGGGCATGA :0.821835750551094
                                                ^
T       G       G       T       C       T       G       C       A       G       A       G       G       G       C       A       T       G       A
2.650   2.650   2.650   1.520   2.650   1.450   1.970   2.850   0.292   0.950   1.380   2.850   0.940   2.850   2.850   1.670   2.850   3.020   3.110
                                ^
A  [ 87 167 281  56   8 744  40 107 851   5 333  54  12  56 104 372  82 117 402 ]
C  [291 145  49 800 903  13 528 433  11   0   3  12   0   8 733  13 482 322 181 ]
G  [ 76 414 449  21   0  65 334  48  32 903 566 504 890 775   5 507 307  73 266 ]
T  [459 187 134  36   2  91  11 324  18   3   9 341   8  71  67  17  37 396  59 ]
                  ^
                                  ^
