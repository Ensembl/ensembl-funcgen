#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


#$* is list of bed filepaths

PASS=$1
shift


$EFG_SRC/scripts/import/parse_and_import.pl\
	-name  rollback_test\
	-format SEQUENCING\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-experimental_set rollback_test\
	-port 3306\
	-host ens-genomics1\
	-dbname expset_homo_sapiens_funcgen_49_36k\
	-group efg\
	-cell_type CD4\
	-feature_type H4R3me2\
	-feature_analysis Parzen\
	-parser Bed\
	-vendor SOLEXA\
	-assembly 36\
	-exp_date 2007-07-17\
	-tee\
	-verbose\
	-recover\
	-pass $PASS\
	$*
