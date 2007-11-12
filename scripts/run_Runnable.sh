#!/bin/sh

if [ $# -lt 5 ]; then

	echo "Usage: $0 <host> <port> <user> <password> <runnable>"
	exit;

fi

HOST=$1
shift
PORT=$1
shift
USER=$1
shift
PASS=$1
shift
RUNNABLE=$1
shift
ARGS=$*

# set environment variables
if [ -e $EFG_SRC/scripts/.efg_${RUNNABLE}_Runnable ]; then
	. $EFG_SRC/scripts/.efg_${RUNNABLE}_Runnable
else
	echo "ERROR: Need to provide file named $EFG_SRC/scripts/.efg_${RUNNABLE}_Runnable"
	echo "  containing environment variables for configuration."
	exit;
fi

OUTDIR=$ANALYSIS_WORK_DIR

#bsub -o "$OUTDIR/${RUNNABLE}_runnable.log" -J ${RUNNABLE}_[1-25] \
    $EFG_SRC/scripts/run_Runnable.pl \
    -host $HOST \
    -port $PORT \
    -user $USER \
    -pass $PASS \
    -dbname sg_homo_sapiens_funcgen_48_36j \
    -species homo_sapiens \
    -data_version 46_36h \
    -logic_name $LOGIC_NAME \
    -module $MODULE \
    $ARGS

#    -input_name 1:1,1000000 \
#    -write \

