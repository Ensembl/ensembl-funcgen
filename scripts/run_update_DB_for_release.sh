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
then echo "Must provide at a user and password argument"; echo "Write run script usage in .efg"; exit; fi


#| bos_taurus_funcgen_56_4e                | DONE
#| caenorhabditis_elegans_funcgen_56_200   | DONE
#| canis_familiaris_funcgen_56_2m          | DONE
#| ciona_intestinalis_funcgen_56_2m        | DONE 
#| danio_rerio_funcgen_56_8b               | DONE
#| drosophila_melanogaster_funcgen_56_513a | DONE
#| gallus_gallus_funcgen_56_2m             | DONE
#| homo_sapiens_funcgen_56_37a             | RUNNING
#| macaca_mulatta_funcgen_56_10l           | DONE
#| mus_musculus_funcgen_56_37i             | RUNNING on ens-genomics1 copied to staging
#| ornithorhynchus_anatinus_funcgen_56_1k  | DONE
#| pan_troglodytes_funcgen_56_21l          | DONE on ens-genomics1 copied to staging
#| rattus_norvegicus_funcgen_56_34x        | RUNNING in screen
#| saccharomyces_cerevisiae_funcgen_56_1j  | DONE
#| xenopus_tropicalis_funcgen_56_41n       | DONE


#Need to do this for all by submitting these jobs to the farm.

#WARNING: This loop only works for standard DB names i.e. no prefixes!

echo "Please check edit the script before running, by adding your dbnames and checking the host parameters"
exit;


for db in sus_scrofa_funcgen_57_9a; do
	latin=$(echo $db | sed 's/_funcgen_.*//')
	data_version=$(echo $db | sed 's/.*funcgen_//')


bsub_cmd="bsub -o $HOME/logs/update_DB_for_release.${latin}.out  -e $HOME/logs/update_DB_for_release.${latin}.err -J update_DB_for_release.${latin} -q basement "


$bsub_cmd time perl -w $EFG_SRC/scripts/update_DB_for_release.pl\
	-species $latin\
	-port 3306 \
	-host ens-genomics1 \
	-user $USER \
	-data_version $data_version\
	-dbname $db\
	-dnadb_host ens-staging2\
	-check_displayable \
	-pass $PASS \
	-tee \

done
