#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import_type.pl\
		-type FeatureType\
       	-name PolII\
        -dbname homo_sapiens_funcgen_46_36h\
		-class POLYMERASE\
		-description 'RNA Polymerase II'\
        -pass $1\

