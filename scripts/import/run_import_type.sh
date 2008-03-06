#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import/import_type.pl\
		-type FeatureType\
       	-name H3K36me3\
        -dbname update_chip_mus_musculus_funcgen_49_37b\
		-description 'TESTING'\
        -pass $1\
		-class HISTONE\
