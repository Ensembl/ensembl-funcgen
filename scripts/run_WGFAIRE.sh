#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")


$EFG_SRC/scripts/parse_and_import.pl\
	-name HeLaS3_FAIRE_WholeGenome\
	-format tiled\
	-location Hinxton\
	-vendor NIMBLEGEN\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-fasta\
	-port 3306\
	-host ens-genomics1\
	-dbname test_homo_sapiens_funcgen_45_36g\
	-array_set\
	-array_name "0701_HG18"\
	-group efg\
	-tee\
	-data_version 43_36e\
	-pass $1\
	-recover


