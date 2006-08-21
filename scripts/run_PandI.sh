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
	-pass $1