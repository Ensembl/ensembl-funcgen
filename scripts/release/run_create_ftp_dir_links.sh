#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

#$* would be list of builds if none specifed will use DEFAULT build for the update DB

#GetOpts here with


USER=$1
shift
#PASS=$1
#shift

#dbname=homo_sapiens_funcgen_67_37
dbname=mus_musculus_funcgen_68_38
dbhost=ens-staging2
dnadb_host=$dbhost
dnadb_user=$USER


if [ ! $USER ]; then
	echo "Must provide a user argument"; exit; 
fi


job_cmd="perl -w $EFG_SRC/scripts/release/create_ftp_dir_links.pl\
	-dbhost $dbhost\
	-dbuser $USER \
	-dbname $dbname\
	-dnadb_host $dnadb_host\
	-dnadb_user $dnadb_user"

#echo -e $job_cmd

$job_cmd

#-pass $PASS






