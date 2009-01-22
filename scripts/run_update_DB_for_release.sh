#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg



#$* would be list of builds if none specifed will use DEFAULT build for the update DB

USER=$1
shift
PASS=$1
shift

if [ ! $PASS ] || [ ! $USER ]
then echo "Must provide at a user and password argument"; echo "Write run script usage in .efg"; exit; fi


perl -w $EFG_SRC/scripts/update_DB_for_release.pl\
	-species mus_musculus\
	-port 3306\
	-host ens-staging\
	-user $USER\
	-data_version 53_37f\
	-dbname mus_musculus_funcgen_53_37f\
	-check_displayable \
	-pass $PASS\
	-skip_meta_coord\
	$@

