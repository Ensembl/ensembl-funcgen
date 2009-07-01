#!/usr/local/bin/bash

#
# AUTHOR:
# Stefan Graf, Ensembl Functional Genomics
#
# LICENCE:
# This code is distributed under an Apache style licence. Please see
# http://www.ensembl.org/info/about/code_licence.html for details.
#
# CONTACT
# Please post comments/questions to the Ensembl development list
# <ensembl-dev@ebi.ac.uk>
#

# Check if environment variables are set
if [ ! $MAQ_WORKDIR ]; then
    echo "ERROR: Need to set and export variable MAQ_WORKDIR to define"
    echo "working directory."
    exit 1
fi

if [ $# -lt 3 ]; then
    echo "ERROR: Arguments missing!"
    echo "   Usage: run_build_profile.sh <fragment length> <bin_size> <file(s)>"
    exit
fi

FRAGLEN=$1
shift

BIN_SIZE=$1
shift

NOF=$#
echo "Submitting job array with $NOF file(s) to the farm"

QUEUE='64bit';

LOG="${f}_build_profile.log";

NAME="build_profile_${BIN_SIZE}";

eval bsub -oo "${NAME}_%I.log" -J "${NAME}[1-${NOF}]" -q "$QUEUE" \
    ${MAQ_BINDIR}/build_profile.pl -bin_size ${BIN_SIZE} -frag_length $FRAGLEN $*
