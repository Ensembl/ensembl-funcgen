#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import/import_type.pl\
	#Mandatory
	-type        "type"\                #e.g. FeatureType CellType Analysis 
    -name        "feature_name"\        #e.g. H3K36me3
    -dbname      "your_db_name"\        #e.g. your_mus_musculus_funcgen_51_37d
	-species     "latin_name"\          #e.g mus_musculus
	-description "wordy description of feature"\
    -pass $PASS\
	-class       "class of FeatureType"\ #e.g. HISTONE
	$ARGS
	#Optional
