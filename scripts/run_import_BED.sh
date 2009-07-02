#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


#$* is list of bed filepaths

PASS=$1
shift


$EFG_SRC/scripts/import/parse_and_import.pl\
	-name  LMI_CD4_H3K4ac\
	-format SEQUENCING\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-group efg\
	-species homo_sapiens\
	-experimental_set LMI_CD4_H3K4ac\
	-port 3306\
	-ucsc_coords\
	-registry_host ens-staging\
	-registry_user ensro\
	-host ens-genomics1\
	-dbname nj_importer_homo_sapiens_funcgen_55_37\
	-cell_type CD4\
	-feature_type H3K4ac\
	-feature_analysis WindowInterval\
	-parser Bed\
	-vendor SOLEXA\
	-assembly 37\
	-tee\
	-recover\
	-pass $PASS\
	$*
