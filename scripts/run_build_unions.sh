#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/build_union_features.pl\
       	-focus_sets Wiggle_H3K4me3,Wiggle_H3K79me3\
		-union_sets Wiggle_H3K27me3,Wiggle_H3K36me3,Wiggle_H3K20me3,Wiggle_H3K79me3,Wiggle_H3K9me3\
        -species homo_sapiens\
        -port 3306\
        -host ens-genomics2\
		-multiplex\
        -dbname union_homo_sapiens_funcgen_45_36g\
		-dump_features\
		-data_version 44_36f\
        -pass $1\
		-out_dir ~/union_test
        $ARGS
