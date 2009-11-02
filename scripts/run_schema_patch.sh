#!/bin/sh


#Why can't this find SRC in the subshell?
#.bashrc variables are lost if we just use sh above, but efg.config vars are not?

#$* would be list of builds if none specifed will use DEFAULT build for the update DB

USER=$1
shift
PASS=$1
shift
SCHEMA_VERSION=$1
shift

usage="e.g.\trun_schema_patch.sh USER PASS SCHEMA_VERSION"

if [ ! $USER ]
then echo -e "Must provide a USER argument[0]\n$usage"; exit; fi

if [ ! $PASS ]
then echo -e "Must provide a PASS argument[1]\n$usage"; exit; fi

if [ ! $SCHEMA_VERSION ]
then echo -e "Must provide a SCHEMA_VERSION argument[2]\n$usage"; exit; fi




for host in ens-staging1 ens-staging2; do


#Will this apply manual patch?

perl -w $SRC/ensembl/misc-scripts/schema_patch.pl\
	--host $host \
	--port 3306 \
	--user $USER \
	--pass $PASS \
	--pattern %_funcgen_57_% \
	--patch_funcgen_database \
	--schema 57 \
	--dry_run 1 \
	--interactive 0\
	--logfile ~/logs/schema_patch.${host}_57.log \
	$@
done
