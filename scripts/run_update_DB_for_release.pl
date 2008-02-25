#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg



#$* would be list of builds if none specifed will use DEFAULT build for the update DB

PASS=$1
shift

if [ ! $PASS ]
then echo "Must provide at least a password argument"; echo "Write run script usage in .efg"; exit; fi




$EFG_SRC/scripts/update_DB_for_release.pl\
  -species homo_sapiens\
  -port 3306\
  -host ens-genomics1\
  -user ensadmin\
  -data_version 49_36k\
  -dbname expset_homo_sapiens_funcgen_49_36k\
  -pass $PASS\
  $*
