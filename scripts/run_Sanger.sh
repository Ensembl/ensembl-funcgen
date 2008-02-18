#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/parse_and_import.pl\
	-name  H3K4me1-GM06990-new\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-fasta\
	-port 3306\
	-host ens-genomics1\
	-dbname sanger_test_43_36e\
	-array_name "ENCODE3.1.1"\
	-group efg\
	-vendor sanger\
	-assembly 36\
	-exp_date 2006-11-02\
	-tee\
	-verbose\
	-pass $1\
	-recover\
	$ARGS
