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

#echo "Please check edit the script before running, by adding your dbnames and checking the host parameters"
#exit;


#host='ens-genomics2'
#host='ens-staging1'
#dnadb_host='ens-staging1'

#host='ens-genomics2'
dnadb_host='ens-staging2'

host=$dnadb_host


#dbs="homo_sapiens_funcgen_72_37 bos_taurus_funcgen_72_31  caenorhabditis_elegans_funcgen_72_235  canis_familiaris_funcgen_72_31    ciona_intestinalis_funcgen_72_3   danio_rerio_funcgen_72_9  drosophila_melanogaster_funcgen_72_546  gallus_gallus_funcgen_72_4  macaca_mulatta_funcgen_72_10"
#dbs="homo_sapiens_funcgen_72_37"

#dbs="mus_musculus_funcgen_72_38"
#dbs="mus_musculus_funcgen_72_38 ornithorhynchus_anatinus_funcgen_72_1 oryctolagus_cuniculus_funcgen_72_3 pan_troglodytes_funcgen_72_214 rattus_norvegicus_funcgen_72_5 saccharomyces_cerevisiae_funcgen_72_4 sus_scrofa_funcgen_72_102 xenopus_tropicalis_funcgen_72_42"



#### SKIP EG DBs ??????? #############


echo "dbs $dbs"
dnadb_port=3306
dnadb_user=$USER
dnadb_pass=$PASS

for db in $dbs; do
	echo -e "\n\n::\tUpdating ${host}:${db}"
	#Put this in the log
	latin=$(echo $db | sed 's/_funcgen_.*//')
	latin=$(echo $latin | sed 's/dev_//')
	data_version=$(echo $db | sed 's/.*funcgen_//')

	bsub_cmd="bsub -o $HOME/logs/update_DB_for_release.${latin}_${data_version}.out  -e $HOME/logs/update_DB_for_release.${latin}_${data_version}.err -J update_DB_for_release.${latin}_${data_version} -q long -R\"select[mem>2000] rusage[mem=2000]\" -M 2000000"


	job_cmd="perl -w $EFG_SRC/scripts/release/update_DB_for_release.pl\
	-species $latin\
	-port 3306 \
	-host $host\
	-user $USER \
	-data_version $data_version\
	-dbname $db\
	-dnadb_host $dnadb_host\
	-dnadb_user $dnadb_user\
	-dnadb_pass $dnadb_pass\
    -dnadb_port $dnadb_port\
	-check_displayable \
	-pass $PASS -no_log $*"

	echo -e "$bsub_cmd $job_cmd"

	#omit -no_log if running locally
	#echo through bash to avoid -R anomalies VOODOO!

	#$job_cmd

	echo "$bsub_cmd $job_cmd" | bash
done

#-skip_meta_coord


exit;





