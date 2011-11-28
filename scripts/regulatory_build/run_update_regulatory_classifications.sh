#!/usr/local/bin/bash

USAGE="Usage: $0 <password> [ script params ]"

if [ $# -lt 1 ]; then
	echo $USAGE
	exit;
fi

PASS=$1
shift
ARGS=$*
 

USER=ensadmin
PORT=3306
SPECIES=mus_musculus
DATA_VERSION=65_37
DB_NAME="dev_${SPECIES}_funcgen_${DATA_VERSION}"
DB_PREFIX="annotation_dev_${SPECIES}_funcgen_${DATA_VERSION}"
HOST=ens-genomics1
OUT_DIR="/lustre/scratch103/ensembl/funcgen/output/${DB_NAME}/"

if [[ ! -d $OUT_DIR ]]; then
	mkdir $OUT_DIR;
fi

LOG_FILE="${OUT_DIR}/update_regulatory_classifications.$$.log"

#if [[ $SPECIES = homo_* ]] || [[ $SPECIES = human ]]; then
#	SPECIES=homo_sapiens
#	HOST='ens-genomics2'
#elif [[ $SPECIES = mus_* ]] || [[ $SPECIES = mouse ]]; then
#	SPECIES=mus_musculus
#	HOST='ens-genomics1'
#else
#	echo -e "$SPECIES not recognised, unable to set default host\n$USAGE"
#	exit
#fi

#echo "Species recognised as: $SPECIES"
#echo "Setting DB host to:    $HOST"


job_cmd="perl -w $EFG_SRC/scripts/regulatory_build/update_regulatory_classifications.pl\
	-dbname $DB_NAME \
 	-host   $HOST \
 	-user   $USER \
   	-pass   $PASS \
    -species $SPECIES\
    -dbprefix $DB_PREFIX\
    -logfile $LOG_FILE\
  $*
"
echo $job_cmd
$job_cmd
  
#  -report_only\
