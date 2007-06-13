#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import_type.pl\
		-type FeatureType\
       	-name H3ac\
        -dbname nath_MA2C_3042_mus_musculus_funcgen_46_36g\
		-class HISTONE\
		-description 'Histone 3 Acetylation'\
        -pass $1\

