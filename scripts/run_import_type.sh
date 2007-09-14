#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import_type.pl\
		-type FeatureType\
       	-name H4R3me2\
        -dbname expset_homo_sapiens_funcgen_47_36i\
		-class HISTONE\
		-description 'TESTING'\
        -pass $1\

