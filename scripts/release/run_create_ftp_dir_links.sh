#!/bin/sh

USER=$1
shift



dbname=your_binomial_speciesname_funcgen_SCHEMA_BUILD
dbhost=your_funcgen_mysql_host

#Optional
dnadb_host=your_dnadb_mysql_host
dnadb_user=$USER

nfs_root='/nfs/root/staging/source/data/directory/species/assmebly';
ftp_root='/ftp/mirror/root/directory'

if [ ! $USER ]; then
	echo "Must provide a user argument"; exit; 
fi


job_cmd="perl -w $EFG_SRC/scripts/release/create_ftp_dir_links.pl\
  -nfs_root $nfs_root\
  -ftp_root $ftp_root\
	-dbhost $dbhost\
	-dbuser $USER \
	-dbname $dbname\
	-dnadb_host $dnadb_host\
	-dnadb_user $dnadb_user $@"

$job_cmd







