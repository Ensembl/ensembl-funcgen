#!/bin/sh

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import_array_from_fasta.pl\
		-format 'CUSTOM'\
		-vendor 'NIMBLEGEN'\
		-array 'MA2C_3042'\
		-dbname nath_MA2C_3042_mus_musculus_funcgen_46_36g\
        -pass $1\
        -file /lustre/work1/ensembl/nj1/affy/nath_MA2C_3042_mus_musculus_funcgen_46_36g/MA2C_3042_probe.fasta\
		-outdir /lustre/work1/ensembl/nj1/affy/nath_MA2C_3042_mus_musculus_funcgen_46_36g\
		$ARGS
