#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

$EFG_SRC/scripts/parse_and_import.pl\
	-name  H3K4me3-GM06990\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-fasta\
	-port 3306\
	-host ens-genomics1\
	-dbname homo_sapiens_funcgen_43\
	-array_name "ENCODE3.1.1"\
	-group efg\
	-vendor sanger\
	-data_version 41_36c\
	-exp_date 2006-11-02\
	-verbose 2\
	-pass $1\
	-recover


