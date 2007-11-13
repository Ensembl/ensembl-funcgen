#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

PASS=$1
shift




perl $EFG_SRC/scripts/external_features/load_external_features.pl\
	-type  cisred\
	-file  /nfs/acari/nj1/src/ensembl-functgenomics/scripts/external_features/cisRED/cisred_9_1_motifs4Ensembl.txt\
	-species homo_sapiens\
	-port 3306\
	-user ensadmin\
	-host ens-genomics1\
	-clobber 1\
	-dbname  external_homo_sapiens_funcgen_48_36j\
	-pass $PASS
