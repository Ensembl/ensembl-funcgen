#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

$EFG_SRC/scripts/parse_and_import.pl \
	-name Nimblegen_CHIP2_data\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-fasta\
    -group efg\
	-data_version 40_36b\
	-verbose 2\
	-pass $1