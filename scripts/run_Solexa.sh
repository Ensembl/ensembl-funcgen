#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


#$* is list of bed filepaths

PASS=$1
shift


$EFG_SRC/scripts/import/parse_and_import.pl\
	-name  NP_H3K4me3\
	-format SEQUENCING\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-group efg\
	-species mus_musculus\
	-experimental_set Mikkelsen_NPC_H3K4me3\
	-port 3306\
	-is_ucsc\
	-host ens-genomics1\
	-dbname mus_musculus_funcgen_50_37c\
	-cell_type NPC\
	-feature_type H3K4me3\
	-feature_analysis WindowInterval\
	-parser Bed\
	-vendor SOLEXA\
	-assembly 37\
	-tee\
	-recover\
	-pass $PASS\
	$*
