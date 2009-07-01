#!/usr/local/bin/bash

if [ $# -lt 3 ]; then
    echo "USAGE: $0 <queue> <prefix> <file(s)>"
    exit 1
fi

QUEUE=$1
shift

PREFIX=$1
shift

NAME="load_bed"
JOB="${NAME}[1-$#]%4"
LOG="${NAME}.log"
RES="rusage[mem=2000]"
#BSUB_RESOURCES='"select[type==X86_64 && mem>3000] rusage[mem=3000]"'
#BSUB_RESOURCES='"select[type==X86_64 && mem>6000] rusage[mem=6000]"'


#cat <<EOF

eval bsub -oo "$LOG" -J "$JOB" -q "$QUEUE" -R "$RES" \
    ${EFG_SRC}/scripts/DAS/load_bed.pl -prefix $PREFIX $@

#EOF

