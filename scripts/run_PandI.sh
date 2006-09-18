#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

$EFG_SRC/scripts/parse_and_import.pl\
	-name Stunnenberg_all_OID_1963\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-fasta\
	-port 3362\
	-host ecs2\
	-dbname test_homo_sapiens_funcgen_41_36b\
	-array_set\
	-array_name "2005-05-10_HG17Tiling_Set"\
	-group efg\
	-data_version 40_36b\
	-verbose 2\
	-pass $1\
	-recover
