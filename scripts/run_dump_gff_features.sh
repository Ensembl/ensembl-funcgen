#!/bin/sh

PASS=$1
shift




$EFG_SRC/scripts/export/dump_gff_features.pl\
        -species homo_sapiens\
        -port 3306\
        -dbhost ens-genomics1\
        -dbname homo_sapiens_funcgen_50_36l\
		-cdbname homo_sapiens_core_50_36l\
        -feature_set 'RegulatoryFeatures' \
		-outdir $HOME/data/efg/reg_build/v50 \
		-pass $1		
