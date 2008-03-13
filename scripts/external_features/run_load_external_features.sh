#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

#file for cisred is always motif file, search regions will be guessed from filename

perl $EFG_SRC/scripts/external_features/load_external_features.pl\
	-type  redfly\
	-species drosophila_melanogaster\
	-port 3306\
	-user ensadmin\
	-host ens-genomics1\
	-clobber 1\
	-dbname not_patched_drosophila_melanogaster_funcgen_50_54a\
	-cdbname drosophila_melanogaster_core_49_54\
	-pass $1



#	-file  /nfs/acari/nj1/src/ensembl-functgenomics/scripts/external_features/VISTA/vista.txt\
