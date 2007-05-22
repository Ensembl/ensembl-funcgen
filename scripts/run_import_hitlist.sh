#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/convert_file_to_features.pl\
       	-feature_type H3ac\
		-format HitList\
        -species homo_sapiens\
        -port 3306\
        -host ens-genomics2\
        -dbname new_homo_sapiens_funcgen_45_36g\
        -set_name "H3ac_claes"\
        -data_version 25_34e\
		-is_ucsc 1\
		-clobber\
        -pass $1\
        -file /lustre/scratch1/ensembl/prf1/h3ac_uu/H3ac.track\
		$ARGS
