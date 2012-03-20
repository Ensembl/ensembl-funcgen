#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#1. DataSets
#Datasets are a meta container for data, grouping the data obtained by an analysis (FeatureSets) to the underlying raw data (ResultSets).
#Create a script which fetches all available DataSets for Human. How many are there?
#Now get the 'RegulatoryFeatures:MultiCell' data set and print the display label of the product feature set and all the supporting sets.

#Grab the eFG adaptor
my $dsa = $registry->get_adaptor('Human', 'funcgen', 'dataset');

my @human_data_sets = @{$dsa->fetch_all()};
print "There are ".scalar(@human_data_sets)." datasets for Human\n";

#foreach my $dataset (@human_data_sets){
#	print $dataset->name."\t".$dataset->cell_type->name."\t".$dataset->feature_type->name."\n";
#}


my $dset = $dsa->fetch_by_name('RegulatoryFeatures:MultiCell');


#Print the product FeatureSet name
print "Product FeatureSet:\t".$dset->product_FeatureSet->name."\n";

#Now print all the supporting sets
my @ssets = @{$dset->get_supporting_sets};
print "Found ".scalar(@ssets)." supporting sets\n";

foreach my $sset(@ssets){
  print "Supporting ".ucfirst($sset->feature_class)."Feature set:\t".$sset->display_label."\n";

}

__END__

>perl sets_1.pl

There are 750 datasets for Human
Product FeatureSet:     RegulatoryFeatures:MultiCell
Found 245 supporting sets
Supporting AnnotatedFeature set:        DNase1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - NHEK Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - GM06990 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - NHEK Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - K562 Enriched Sites
Supporting AnnotatedFeature set:        FAIRE - NHEK Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - IMR90 Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HMEC Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - HSMM Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - NH-A Enriched Sites
Supporting AnnotatedFeature set:        DNase1 - IMR90 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - CD4 Enriched Sites
Supporting AnnotatedFeature set:        Gabp - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cjun - K562 Enriched Sites
Supporting AnnotatedFeature set:        Jund - K562 Enriched Sites
Supporting AnnotatedFeature set:        Max - K562 Enriched Sites
Supporting AnnotatedFeature set:        Nfe2 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Srf - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - K562 Enriched Sites
Supporting AnnotatedFeature set:        Nrsf - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cfos - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - NHEK Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Gabp - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Srf - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Nrsf - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Cfos - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Cfos - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Jund - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Max - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Max - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - GM06990 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        CTCF - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - NHEK Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        BCL3 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        EBF - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Egr1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        IRF4 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        p300 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Pax5 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Sin3Ak20 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        USF1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        ZBTB33 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Nrsf - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        TAF1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        FOSL2 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        HEY1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        RXRA - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        ZBTB33 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Egr1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        HEY1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        PU1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Sin3Ak20 - K562 Enriched Sites
Supporting AnnotatedFeature set:        SIX5 - K562 Enriched Sites
Supporting AnnotatedFeature set:        TAF1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        USF1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Yy1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Ap2alpha - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Ap2gamma - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        E2F1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        E2F4 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        E2F6 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Tr4 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Bdp1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Brf1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Rad21 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Sirt6 - K562 Enriched Sites
Supporting AnnotatedFeature set:        XRCC4 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        CTCF - NHEK Enriched Sites
Supporting AnnotatedFeature set:        BATF - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        BCL11A - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        NFKB - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Rad21 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Tr4 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        ZZZ3 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        BAF155 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        BAF170 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Bdp1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Brf1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Brf2 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Brg1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Cjun - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Ini1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Jund - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Nrf1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        RPC155 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        TFIIIC-110 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        Cjun - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        Max - HUVEC Enriched Sites
Supporting AnnotatedFeature set:        ATF3 - K562 Enriched Sites
Supporting AnnotatedFeature set:        E2F4 - K562b Enriched Sites
Supporting AnnotatedFeature set:        E2F6 - K562b Enriched Sites
Supporting AnnotatedFeature set:        Gata1 - K562b Enriched Sites
Supporting AnnotatedFeature set:        Gata2 - K562b Enriched Sites
Supporting AnnotatedFeature set:        Brf2 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Brg1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        SETDB1 - K562b Enriched Sites
Supporting AnnotatedFeature set:        Tr4 - K562b Enriched Sites
Supporting AnnotatedFeature set:        Yy1 - K562b Enriched Sites
Supporting AnnotatedFeature set:        Znf263 - K562b Enriched Sites
Supporting AnnotatedFeature set:        ZNF274 - K562b Enriched Sites
Supporting AnnotatedFeature set:        GTF2B - K562 Enriched Sites
Supporting AnnotatedFeature set:        Ini1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        NELFe - K562 Enriched Sites
Supporting AnnotatedFeature set:        TFIIIC-110 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Pbx3 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        POU2F2 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        PU1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        SP1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        TAF1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Tcf12 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Gabp - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        TAF1 - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        BHLHE40 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Gabp - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Jund - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        p300 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Sin3Ak20 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        USF1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        SP1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        CBP - CD4 Enriched Sites
Supporting AnnotatedFeature set:        HDAC2 - CD4 Enriched Sites
Supporting AnnotatedFeature set:        HDAC3 - CD4 Enriched Sites
Supporting AnnotatedFeature set:        HDAC6 - CD4 Enriched Sites
Supporting AnnotatedFeature set:        MOF - CD4 Enriched Sites
Supporting AnnotatedFeature set:        p300 - CD4 Enriched Sites
Supporting AnnotatedFeature set:        PCAF - CD4 Enriched Sites
Supporting AnnotatedFeature set:        Tip60 - CD4 Enriched Sites
Supporting AnnotatedFeature set:        HDAC1 - CD4 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HMEC Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HMEC Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HSMM Enriched Sites
Supporting AnnotatedFeature set:        CTCF - NH-A Enriched Sites
Supporting AnnotatedFeature set:        Cjun - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Cmyc - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Gata2 - K562 Enriched Sites
Supporting AnnotatedFeature set:        HDAC8 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Junb - K562 Enriched Sites
Supporting AnnotatedFeature set:        Jund - K562 Enriched Sites
Supporting AnnotatedFeature set:        NR4A1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        ATF3 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        BCLAF1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        ELF1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        ETS1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        MEF2A - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        MEF2C - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Rad21 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        RXRA - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        SIX5 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        Yy1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        ZEB1 - GM12878 Enriched Sites
Supporting AnnotatedFeature set:        ATF3 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        BCL11A - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        CTCF - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Egr1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        FOSL1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Gabp - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        HDAC2 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Jund - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Nanog - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        p300 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        POU5F1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Rad21 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        RXRA - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Sin3Ak20 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        SIX5 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        SP1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Srf - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        TAF7 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Tcf12 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        USF1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        Yy1 - H1ESC Enriched Sites
Supporting AnnotatedFeature set:        ATF3 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        ELF1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        FOXA1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        FOXA2 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        HDAC2 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        HNF4A - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        HNF4G - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Nrsf - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Rad21 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        SP1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Srf - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        TAF1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Tcf12 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        Yy1 - HepG2 Enriched Sites
Supporting AnnotatedFeature set:        ATF3 - K562 Enriched Sites
Supporting AnnotatedFeature set:        BCL3 - K562 Enriched Sites
Supporting AnnotatedFeature set:        BCLAF1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        CTCFL - K562 Enriched Sites
Supporting AnnotatedFeature set:        E2F6 - K562 Enriched Sites
Supporting AnnotatedFeature set:        ELF1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        ETS1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        FOSL1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Gata2 - K562 Enriched Sites
Supporting AnnotatedFeature set:        HDAC2 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Max - K562 Enriched Sites
Supporting AnnotatedFeature set:        MEF2A - K562 Enriched Sites
Supporting AnnotatedFeature set:        Rad21 - K562 Enriched Sites
Supporting AnnotatedFeature set:        SP2 - K562 Enriched Sites
Supporting AnnotatedFeature set:        TAF7 - K562 Enriched Sites
Supporting AnnotatedFeature set:        THAP1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        Yy1 - K562 Enriched Sites
Supporting AnnotatedFeature set:        ZBTB33 - K562 Enriched Sites
Supporting AnnotatedFeature set:        ZBTB7A - K562 Enriched Sites
Supporting AnnotatedFeature set:        Nrsf - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - HeLa-S3 Enriched Sites
Supporting AnnotatedFeature set:        CTCF - GM12878 Enriched Sites
