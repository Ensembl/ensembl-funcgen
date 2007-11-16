#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

#file for cisred is always motif file, search regions will be guessed from filename

perl $EFG_SRC/scripts/external_features/load_external_features.pl\
	-type  vista\
	-file  /nfs/acari/nj1/src/ensembl-functgenomics/scripts/external_features/VISTA/vista.txt\
	-species homo_sapiens\
	-port 3306\
	-user ensadmin\
	-host ens-genomics1\
	-clobber 1\
	-dbname homo_sapiens_funcgen_48_36j\
	-pass $1
