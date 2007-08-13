#!/bin/sh

PASS=$1
shift
ARGS=$*

#bsub -o "/lustre/scratch1/ensembl/graef/RegBuild/v46/RegulatoryBuild.log" -J RegBuild \

$EFG_SRC/scripts/build_regulatory_features.pl-2007-08-02 \
    -overlap_sets Overlap_GM06990_DNASE_IMPORT::CD4_CTCF:CD4_H2AZ:CD4_H2BK5me1:CD4_H3K27me1:CD4_H3K27me2:CD4_H3K27me3:CD4_H3K36me1:CD4_H3K36me3:CD4_H3K4me1:CD4_H3K4me2:CD4_H3K4me3:CD4_H3K79me1:CD4_H3K79me2:CD4_H3K79me3:CD4_H3K9me1:CD4_H3K9me2:CD4_H3K9me3:CD4_H3R2me1:CD4_H3R2me2:CD4_H4K20me1:CD4_H4K20me3:CD4_H4R3me2:CD4_PolII:GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3:Wiggle_H4K20me3,Overlap_CD4_DNASE_IMPORT::CD4_CTCF:CD4_DNASE_IMPORT:CD4_H2AZ:CD4_H2BK5me1:CD4_H3K27me1:CD4_H3K27me2:CD4_H3K27me3:CD4_H3K36me1:CD4_H3K36me3:CD4_H3K4me1:CD4_H3K4me2:CD4_H3K4me3:CD4_H3K79me1:CD4_H3K79me2:CD4_H3K79me3:CD4_H3K9me1:CD4_H3K9me2:CD4_H3K9me3:CD4_H3R2me1:CD4_H3R2me2:CD4_H4K20me1:CD4_H4K20me3:CD4_H4R3me2:CD4_PolII:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3:Wiggle_H4K20me3\
	-default_set Overlap_CD4_DNASE_IMPORT::CD4_CTCF:CD4_DNASE_IMPORT:CD4_H2AZ:CD4_H2BK5me1:CD4_H3K27me1:CD4_H3K27me2:CD4_H3K27me3:CD4_H3K36me1:CD4_H3K36me3:CD4_H3K4me1:CD4_H3K4me2:CD4_H3K4me3:CD4_H3K79me1:CD4_H3K79me2:CD4_H3K79me3:CD4_H3K9me1:CD4_H3K9me2:CD4_H3K9me3:CD4_H3R2me1:CD4_H3R2me2:CD4_H4K20me1:CD4_H4K20me3:CD4_H4R3me2:CD4_PolII:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3:Wiggle_H4K20me3\
    -species homo_sapiens\
    -port 3306\
    -host ens-genomics1\
    -dbname homo_sapiens_funcgen_46_36h\
		-dump_features\
		-data_version 45_36g\
    -pass $PASS\
		-out_dir ~/graef/RegBuild/v46\
    $ARGS


#		-out_dir /lustre/scratch1/ensembl/graef/RegBuild/v46\
#		-clobber\
#		-write_features\

# -overlap_sets Overlap_Nessie_NG_STD_2_ctcf_ren_BR1::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3,Overlap_Wiggle_H3K4me3::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3,Overlap_GM06990_DNASE_IMPORT::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3\
#	-default_set Overlap_GM06990_DNASE_IMPORT::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3\
