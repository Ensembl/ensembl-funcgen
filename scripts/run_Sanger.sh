#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

ARGS=$@

PASS=$1

ARGS=$(echo $ARGS | sed "s/$PASS //")


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
	-data_version 42_36d\
	-exp_date 2006-11-02\
	-verbose 2\
	-pass $1\
	-recover $ARGS


