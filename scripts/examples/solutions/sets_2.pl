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

#2. FeatureSets
#Feature Sets hold processed data or features i.e. peak calls or the output of a high level analysis e.g. the Regulatory Build.
#Print the name of the feature sets for the Human 'GM12878' cell type.
#Print the name of the feature sets for the Human 'CTCF' feature type.
#Is the Human FeatureSet 'VISTA enhancer set' associated to any cell type or feature type?
#Is the VISTA FeatureSet associated to a DataSet?
#Hint: Use fetch_by_product_FeatureSet.
#Any idea why this is so?
#Is the Human FeatureSet 'RegulatoryFeatures:MultiCell' associated to a DataSet?
#Any idea why this is so?

#Grab the eFG adaptor
my $fsa = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my $cta = $registry->get_adaptor('Human', 'funcgen', 'celltype');


#Print GM12878 cell type details
my @GM12878_feature_sets = @{$fsa->fetch_all_by_CellType($cta->fetch_by_name('GM12878'))};
print "There are ".scalar(@GM12878_feature_sets)." feature sets for GM12878:\n";

foreach my $featureset (@GM12878_feature_sets){
	print "\t".$featureset->name."\n";
}


=pod

There are 72 feature sets for GM12878:
        GM12878_CTCF_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_H3K9ac_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_H4K20me1_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_CTCF_ENCODE_Uta_SWEmbl_R015_D150
        GM12878_H3K36me3_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_H3K4me3_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_H3K4me2_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_Gabp_ENCODE_Hudsonalpha_SWEmbl_R015_D150
        GM12878_Srf_ENCODE_Hudsonalpha_SWEmbl_R015_D150
        GM12878_Nrsf_ENCODE_Hudsonalpha_SWEmbl_R015_D150
        GM12878_PolII_ENCODE_Hudsonalpha_SWEmbl_R015_D150
        GM12878_H3K27me3_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_H3K27ac_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_H3K4me1_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_DNase1_ENCODE_Duke_SWEmbl_R0025_D150
        GM12878_Cfos_ENCODE_Yale_SWEmbl_R015_D150
        GM12878_Jund_ENCODE_Yale_SWEmbl_R015_D150
        GM12878_PolII_ENCODE_Yale_SWEmbl_R015_D150
        GM12878_PolIII_ENCODE_Yale_SWEmbl_R015_D150
        GM12878_Max_ENCODE_Yale_SWEmbl_R015_D150
        GM12878_DNase1_ENCODE_Uw_SWEmbl_R0025_D150
        GM12878_H3K36me3_ENCODE_Uw_SWEMBL_R015
        GM12878_BCL3_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_EBF_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_Egr1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_IRF4_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_p300_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_Pax5_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_Sin3Ak20_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_USF1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_ZBTB33_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_FAIRE_ENCODE_UNC_SWEMBL_R0025
        GM12878_Yy1_ENCODE_Yale_SWEMBL_R015
        GM12878_Cmyc_ENCODE_Uta_SWEMBL_R015
        GM12878_PolII_ENCODE_Uta_SWEMBL_R015
        GM12878_H3K27me3_ENCODE_Uw_SWEMBL_R015
        GM12878_BATF_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_BCL11A_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_NFKB_ENCODE_Yale_SWEMBL_R015
        GM12878_Rad21_ENCODE_Yale_SWEMBL_R015
        GM12878_Tr4_ENCODE_Yale_SWEMBL_R015
        GM12878_ZZZ3_ENCODE_Yale_SWEMBL_R015
        GM12878_Pbx3_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_POU2F2_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_PU1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_SP1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_TAF1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_Tcf12_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_H3K36me3_ENCODE_Broad_ccat_histone
        GM12878_H3K36me3_ENCODE_Uw_ccat_histone
        GM12878_Cjun_ENCODE_Yale_SWEMBL_R015
        GM12878_Cmyc_ENCODE_Yale_SWEMBL_R015
        GM12878_H2AZ_ENCODE_Broad_SWEMBL_R015
        GM12878_H3K79me2_ENCODE_Broad_SWEMBL_R015
        GM12878_H3K9me3_ENCODE_Broad_SWEMBL_R015
        GM12878_ATF3_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_BCLAF1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_ELF1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_ETS1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_MEF2A_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_MEF2C_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_Rad21_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_RXRA_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_SIX5_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_Yy1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_ZEB1_ENCODE_Hudsonalpha_SWEMBL_R015
        GM12878_CTCF_ENCODE_Uw_SWEMBL_R015
        GM12878_H3K4me3_ENCODE_Uw_SWEMBL_R015
        GM12878_H3K27me3_ENCODE_Broad_CCAT_HISTONE
        GM12878_H3K27me3_ENCODE_Uw_CCAT_HISTONE
        Segmentation:GM12878
        RegulatoryFeatures:GM12878

=cut




#Print CTCF feature type details
my $fta = $registry->get_adaptor('Human', 'funcgen', 'featuretype');
my @CTCF_feature_sets = @{$fsa->fetch_all_by_FeatureType($fta->fetch_by_name('CTCF'))};
print "\n\nThere are ".scalar(@CTCF_feature_sets)." feature sets for CTCF:\n";

foreach my $featureset (@CTCF_feature_sets){
	print "\t".$featureset->name."\n";
}

=pod

There are 28 feature sets for CTCF:
        K562_CTCF_ENCODE_Broad_SWEmbl_R015_D150
        CD4_CTCF_BarskiZhao_PMID17512414_SWEmbl_R015_D150
        HepG2_CTCF_ENCODE_Uta_SWEmbl_R015_D150
        GM12878_CTCF_ENCODE_Broad_SWEmbl_R015_D150
        GM12878_CTCF_ENCODE_Uta_SWEmbl_R015_D150
        NHEK_CTCF_ENCODE_Broad_SWEmbl_R015_D150
        HeLa-S3_CTCF_ENCODE_Uta_SWEmbl_R015_D150
        H1ESC_CTCF_ENCODE_Broad_SWEmbl_R015_D150
        HUVEC_CTCF_ENCODE_Broad_SWEmbl_R015_D150
        HepG2_CTCF_ENCODE_Uw_SWEmbl_R015_D150
        GM06990_CTCF_ENCODE_Uw_SWEmbl_R015_D150
        K562_CTCF_ENCODE_Uta_SWEmbl_R015_D150
        HeLa-S3_CTCF_ENCODE_Uw_SWEMBL_R015
        HUVEC_CTCF_ENCODE_Uw_SWEMBL_R015
        K562_CTCF_ENCODE_Uw_SWEMBL_R015
        NHEK_CTCF_ENCODE_Uw_SWEMBL_R015
        HepG2_CTCF_ENCODE_Broad_SWEMBL_R015
        HUVEC_CTCF_ENCODE_Uta_SWEMBL_R015
        NHEK_CTCF_ENCODE_Uta_SWEMBL_R015
        HMEC_CTCF_ENCODE_Uw_SWEMBL_R015
        HMEC_CTCF_ENCODE_Broad_SWEMBL_R015
        HSMM_CTCF_ENCODE_Broad_SWEMBL_R015
        NH-A_CTCF_ENCODE_Broad_SWEMBL_R015
        H1ESC_CTCF_ENCODE_Uta_SWEMBL_R015
        H1ESC_CTCF_ENCODE_Hudsonalpha_SWEMBL_R015
        HepG2_CTCF_ENCODE_Hudsonalpha_SWEMBL_R015
        HeLa-S3_CTCF_ENCODE_Broad_SWEMBL_R015
        GM12878_CTCF_ENCODE_Uw_SWEMBL_R015


=cut


#Print VISTA sets details
my $vista_set = $fsa->fetch_by_name('VISTA enhancer set');
print "\n\n".$vista_set->display_label."\n";

my $cell_type = ($vista_set->cell_type) ? $vista_set->cell_type->name : 'NON-DEFINED';
print "\tCell type:\t".$cell_type;

my $ftype = ($vista_set->feature_type) ? $vista_set->feature_type->name : 'NON-DEFINED';
print "\tFeatureType:\t$ftype\n";

my $dsa     = $registry->get_adaptor('Human', 'funcgen', 'dataset');
my $dset    = $dsa->fetch_by_product_FeatureSet($vista_set);
my $dataset = ($dset) ? $dset->name : 'NON-DEFINED';
print "\tDataSet:\t$dataset\n";

=pod


VISTA Enhancers
        Cell type:      NON-DEFINED     FeatureType:    VISTA Target
        DataSet:        NON-DEFINED

=cut


__END__
