#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import_type.pl\
		-type FeatureType\
       	-name H3ac\
        -dbname nath_homo_sapiens_funcgen_46_36h\
		-class HISTONE\
		-description 'Testing'\
        -pass $1\

