#!/usr/bin/perl
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


__END__


>perl features_4.pl

Relative affinity of the MA0139.1 motif TGGTCTGCAGAGGGCATGA :0.863462959414993
Relative affinity of changed motif      TGGTGTGCAGAGGGCATGA :0.747783950261183
                                            ^
T       G       G       T       C       T       G       C       A       G       A       G       G       G       C       A       T       G       A
2.650   2.650   2.650   1.520   2.650   1.450   1.970   2.850   0.292   0.950   1.380   2.850   0.940   2.850   2.850   1.670   2.850   3.020   3.110
                                ^
