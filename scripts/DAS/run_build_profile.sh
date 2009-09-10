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
# GetOpts queue
# Do we need to move the QUEUE params from pipeline.config to efg.config?
# This will inevitably be superceded by matrix implementation?
# MAQ_WORKDIR isn't used?
# Take SAM format, needs import parser? This would be similar import to GFF
# Optionally dump to file or load to new table in separate or eFG DB?

# Check if environment variables are set
#if [ ! $MAQ_WORKDIR ]; then
#    echo "ERROR: Need to set and export variable MAQ_WORKDIR to define"
#    echo "working directory."
#    exit 1
#fi

if [ $# -lt 3 ]; then
    echo "ERROR: Arguments missing!"
    echo "   Usage: run_build_profile.sh <fragment length> <bin_size> <file(s)>"
    exit
fi

FRAGLEN=$1
shift

#default is 150???

BIN_SIZE=$1
shift
#Default is 25???
#Plug this into the collection code to generate multiple bin sizes?

NOF=$#
echo "Submitting job array with $NOF file(s) to the farm"

#QUEUE='64bit';
QUEUE='normal'

LOG="${f}_build_profile.log";

#Change this to include filename
NAME="build_profile_${BIN_SIZE}";


eval bsub -oo "${NAME}_%I.log" -J "${NAME}[1-${NOF}]" -q "$QUEUE" \
    $EFG_SRC/scripts/DAS/load_bed_source.pl  -host ens-genomics1 --user ensadmin --dbname nj1_DAS -bin_size 25 -frag_length 150 --files ES_DNase_le1m_reads.bed.gz --names test_reads -reads --pass ensembl --no_load
