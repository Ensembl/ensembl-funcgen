#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


PASS=$1
shift


perl $EFG_SRC/scripts/external_features/load_external_features.pl\
	-type  cisred\
	-species homo_sapiens\
	-port 3306\
	-user ensadmin\
	-host ens-genomics1\
	-clobber 1\
	-dbname test_homo_sapiens_funcgen_51_36m\
	-dnadb_host ens-staging\
    -dnadb_name homo_sapiens_core_51_36m\
	-dnadb_user ensro\
	-pass $PASS\
	$@


#$@ is files, always motifs then search for cisRED? Or should we do some pattern matching
#separate loading of feature_sets, then we can load just one at a time?


#	-file  /nfs/acari/nj1/src/ensembl-functgenomics/scripts/external_features/VISTA/vista.txt\
