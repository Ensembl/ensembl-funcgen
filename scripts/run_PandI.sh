#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg
#x=50

#while [ $x -lt 51 ]
#do


#echo "Importing run $x"

#ln -s /nfs/acari/nj1/data/efg/input/nimblegen/Stunnenberg_all_OID_1963 /nfs/acari/nj1/data/efg/input/nimblegen/${x}_Stunnenberg_all_OID_1963

PASS=$1


time $EFG_SRC/scripts/import/parse_and_import.pl\
	-name LF2_H3K36me3_NCMLS\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species mus_musculus\
	-vendor NIMBLEGEN\
	-cell_type LF2\
	-feature_type H3K36me3\
	-port 3306\
	-host ens-genomics1\
	-dbname update_chip_mus_musculus_funcgen_49_37b\
	-array_set\
	-array_name "2006-07-17_MM8Tiling_Set"\
	-group efg\
	-assembly 37\
	-tee\
	-pass $PASS\
	-recover

#echo "return val is $?"

#x=`expr $x + 1`

#done
