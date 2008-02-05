#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/convert_file_to_features.pl\
       	-feature_type H3K9me3\
		-format Wiggle\
		-is_ucsc 0\
        -species homo_sapiens\
        -port 3306\
        -host ens-genomics2\
        -dbname new_homo_sapiens_funcgen_45_36g\
        -set_name "Wiggle_H3K9me3"\
        -data_version 44_36f\
        -pass $1\
        -file /lustre/scratch1/ensembl/prf1/histones_bcgsc/K9_ht12.wig\
		$ARGS
