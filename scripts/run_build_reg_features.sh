#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/build_regulatory_features.pl\
       	-overlap_sets Overlap_Nessie_NG_STD_2_ctcf_ren_BR1::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3,Overlap_Wiggle_H3K4me3::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3,Overlap_GM06990_DNASE_IMPORT::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3\
		-default_set Overlap_GM06990_DNASE_IMPORT::GM06990_DNASE_IMPORT:Nessie_NG_STD_2_ctcf_ren_BR1:Wiggle_H3K20me3:Wiggle_H3K27me3:Wiggle_H3K36me3:Wiggle_H3K4me3:Wiggle_H3K79me3:Wiggle_H3K9me3\
        -species homo_sapiens\
		-write_features\
        -port 3306\
        -host ens-genomics2\
        -dbname homo_sapiens_funcgen_45_36g\
		-dump_features\
		-data_version 44_36f\
        -pass $1\
		-out_dir /lustre/scratch1/ensembl/nj1/reg_build\
		-clobber\
	    $ARGS
