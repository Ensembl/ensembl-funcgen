#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg



#$* would be list of builds if none specifed will use DEFAULT build for the update DB

#GetOpts here with

#host
#no_farm
#dbname
#user
#pass
#regex?

USER=$1
shift
PASS=$1
shift



if [ ! $PASS ] || [ ! $USER ]
then echo "Must provide at a user and password argument"; exit; fi

#WARNING: This loop only works for standard DB names (i.e. only dev prefix allowed)

echo "Please check edit the script before running, by adding your dbnames and checking the host parameters"
exit;

host='ens-staging2'
dnadb_host='ens-staging2'
dbs='mus_musculus_funcgen_59_37l'


for db in $dbs; do

	latin=$(echo $db | sed 's/_funcgen_.*//')
	latin=$(echo $latin | sed 's/dev_//')
	data_version=$(echo $db | sed 's/.*funcgen_//')

	bsub_cmd="bsub -o $HOME/logs/update_DB_for_release.${latin}_${data_version}.out  -e $HOME/logs/update_DB_for_release.${latin}_${data_version}.err -J update_DB_for_release.${latin}_${data_version} -q basement "


#	job_cmd="$bsub_cmd time perl -w $EFG_SRC/scripts/update_DB_for_release.pl\
job_cmd="time perl -w $EFG_SRC/scripts/release/update_DB_for_release.pl\
	-species $latin\
	-port 3306 \
	-host $host\
	-user $USER \
	-data_version $data_version\
	-dbname $db\
	-dnadb_host $dnadb_host\
	-check_displayable \
	-pass $PASS $*"

	echo -e "\n$job_cmd"
	
	$job_cmd

done

exit


