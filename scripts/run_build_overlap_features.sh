#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")

for i in Wiggle_H3K4me3 GM06990_DNASE_IMPORT Nessie_NG_STD_2_ctcf_ren_BR1; do

    bsub -J "BuildFeatures_$i_[1-25]" \
         -o /lustre/scratch1/ensembl/nj1/reg_build/${i}_%I.log \
		$EFG_SRC/scripts/build_overlap_features.pl \
		-f  $i \
		-t Wiggle_H3K4me3,Wiggle_H3K9me3,Wiggle_H3K20me3,Wiggle_H3K27me3,Wiggle_H3K36me3,Wiggle_H3K79me3,GM06990_DNASE_IMPORT,Nessie_NG_STD_2_ctcf_ren_BR1 \
		-h ens-genomics2 \
		-o $EFG_DATA \
		-p $1 \
		-v 44_36f \
		-d sg_reg_homo_sapiens_funcgen_45_36g -w -c
done

