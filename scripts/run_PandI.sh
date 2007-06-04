#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg
#x=50

#while [ $x -lt 51 ]
#do


#echo "Importing run $x"

#ln -s /nfs/acari/nj1/data/efg/input/nimblegen/Stunnenberg_all_OID_1963 /nfs/acari/nj1/data/efg/input/nimblegen/${x}_Stunnenberg_all_OID_1963

PASS=$1


$EFG_SRC/scripts/parse_and_import.pl\
	-name Stunnenberg_all_OID_1963\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-vendor NIMBLEGEN\
	-fasta\
	-port 3306\
	-host ens-genomics1\
	-dbname maptest_human_funcgen_45_36g\
	-array_set\
	-array_name "2005-05-10_HG17Tiling_Set"\
	-group efg\
	-data_version 44_36f\
	-farm 0\
	-verbose\
	-tee\
	-pass $PASS\
	-recover

#echo "return val is $?"

#x=`expr $x + 1`

#done
