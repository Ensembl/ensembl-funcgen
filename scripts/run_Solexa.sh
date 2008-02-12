#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


#$ARGS is list of bed filepaths

PASS=$1
shift



$EFG_SRC/scripts/import/parse_and_import.pl\
	-name  CD4_H4R3me2\
	-format SEQUENCING\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-experimental_set testing\
	-port 3306\
	-host ens-genomics1\
	-dbname expset_homo_sapiens_funcgen_47_36i\
	-group efg\
	-cell_type CD4\
	-feature_type H4R3me2\
	-vendor SOLEXA\
	-parser Bed\
	-data_version 46_36h\
	-exp_date 2007-07-17\
	-tee\
	-verbose\
	-pass $1\
	-recover\
	$ARGS
