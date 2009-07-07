#!/usr/local/bin/bash

#LICENSE
#
#  Copyright (c) 1999-2009 The European Bioinformatics Institute and
#  Genome Research Limited.  All rights reserved.
#
#  This software is distributed under a modified Apache license.
#  For license details, please see
#
#    http://www.ensembl.org/info/about/code_licence.html
#
#CONTACT
#
#  Please email comments or questions to the public Ensembl
#  developers list at <ensembl-dev@ebi.ac.uk>.
#
#  Questions may also be sent to the Ensembl help desk at
#  <helpdesk@ensembl.org>.

# To do
# Integrate this into efg.env as a function? BuildBedProfile?
# GetOpts
# Do we need to move the QUEUE params from pipeline.config to efg.config?
# This will inevitably be superceded by matrix implementation?
# MAQ_WORKDIR isn't used?

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
