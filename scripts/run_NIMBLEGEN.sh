#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


PASS=$1

$EFG_SRC/scripts/parse_and_import.pl\
	-name "DVD/EXPERIMENT_NAME"\
	-format tiled\
	-location Hinxton\
	-contact "your_email@address.com\
	-species homo_sapiens\
	-fasta\
	-port 3306\
	-host dbhost\
	-dbname "dbname_funcgen_VERSION_BUILD"\
	-array_set\
	-array_name "DESIGN_NAME"\
	-group efg\
	-data_version 41_36c\
	-verbose\
	-tee\
	-pass $PASS\
	-recover

