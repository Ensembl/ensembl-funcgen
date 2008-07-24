#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg
#x=50

#while [ $x -lt 51 ]
#do


#echo "Importing run $x"

#ln -s /nfs/acari/nj1/data/efg/input/nimblegen/Stunnenberg_all_OID_1963 /nfs/acari/nj1/data/efg/input/nimblegen/${x}_Stunnenberg_all_OID_1963

PASS=$1


time $EFG_SRC/scripts/import/parse_and_import.pl\
	-name Stunnenberg_all_OID_1963\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-vendor NIMBLEGEN\
	-cell_type Hela\
	-feature_type H3ac\
	-port 3306\
	-host ens-genomics1\
	-dbname cachetest_homo_sapiens_funcgen_50_36l\
	-array_set\
	-array_name "2005-05-10_HG17Tiling_Set"\
	-group efg\
	-assembly 36\
	-tee\
	-pass $PASS\
	-recover

#echo "return val is $?"

#x=`expr $x + 1`

#done
