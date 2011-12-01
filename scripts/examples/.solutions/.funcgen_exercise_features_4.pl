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

my $rf = $rfa->fetch_by_stable_id('ENSR00001227187');
my $mf = $rf->regulatory_attributes('motif')->[0];

my $slice = $mf->feature_Slice;

my $seq = $mf->seq;

my $bm = $mf->binding_matrix;
print "Relative affinity of the reference motif $seq :".$bm->relative_affinity($seq)."\n";
print "Relative affinity of changed motif CTGCCTACAGAGGTAGCAC :".$bm->relative_affinity("CTGCCTACAGAGGTAGCAC")."\n";

my $cs_adaptor = $registry->get_adaptor("Multi", 'compara', 'ConservationScore');
my $mlss_adaptor = $registry->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
#get method_link_species_set object for gerp conservation scores for mammals
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");

my $display_size = $slice->end - $slice->start + 1; 
my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice, $display_size);
print join("\t",split('',$seq))."\n";
map { print "\t".$_->diff_score; } @$scores; print "\n";
