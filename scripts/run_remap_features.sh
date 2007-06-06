#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/remap_features.pl\
        -species homo_sapiens\
        -port 3306\
		-old_assembly NCBI35\
		-new_assembly NCBI36\
        -host ens-genomics2\
        -dbname new_homo_sapiens_funcgen_45_36g\
        -set_name GM06990_DNASE_IMPORT\
        -new_data_version 44_36f\
		-old_data_version 36_35i\
        -pass $1		
