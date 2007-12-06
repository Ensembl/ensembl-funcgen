#!/bin/sh

PASS=$1
shift


$EFG_SRC/scripts/import_type.pl\
		-type CellType\
       	-name U2OS\
        -dbname tukey_homo_sapiens_funcgen_48_36j\
		-description 'TESTING'\
        -pass $1
#		-class HISTONE\
