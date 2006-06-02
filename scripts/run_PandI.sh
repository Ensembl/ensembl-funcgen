#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

$EFG_SRC/scripts/parse_and_import.pl \
	-instance Nimblegen_CHIP2_data \
	-format tiled\
	-fasta\
    -group efg\
	-pass ensembl