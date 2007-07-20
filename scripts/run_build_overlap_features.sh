#!/bin/sh
USER=$1
PASS=$2

focus="CD4_DNASE_IMPORT GM06990_DNASE_IMPORT"
target="\
Nessie_NG_STD_2_ctcf_ren_BR1,\
CD4_CTCF,\
CD4_H2AZ,\
CD4_H2BK5me1,\
CD4_H3K27me1,\
CD4_H3K27me2,\
CD4_H3K27me3,\
CD4_H3K36me1,\
CD4_H3K36me3,\
CD4_H3K4me1,\
CD4_H3K4me2,\
CD4_H3K4me3,\
CD4_H3K79me1,\
CD4_H3K79me2,\
CD4_H3K79me3,\
CD4_H3K9me1,\
CD4_H3K9me2,\
CD4_H3K9me3,\
CD4_H3R2me1,\
CD4_H3R2me2,\
CD4_H4K20me1,\
CD4_H4K20me3,\
CD4_H4R3me2,\
CD4_PolII,\
Wiggle_H3K27me3,\
Wiggle_H3K36me3,\
Wiggle_H3K4me3,\
Wiggle_H3K79me3,\
Wiggle_H3K9me3,\
Wiggle_H4K20me3"

for f in $focus; do

    bsub -J "BuildFeatures_$f_[1-25]" \
         -o $EFG_DATA/overlap_features/v46/${f}_%I.log \
		$EFG_SRC/scripts/build_overlap_features.pl \
		-f $f \
		-t $target \
		-h ens-genomics1 \
        -u $USER \
        -p $PASS \
		-o $EFG_DATA/overlap_features/v46 \
		-v 45_36g \
		-d homo_sapiens_funcgen_46_36h # -w -c
    
done



