#!/bin/sh

if [ $# -lt 4 ]; then

	echo "Usage: $0 <host> <port> <user> <password>"
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
ARGS=$*

FOCUS=\
CD4_CTCF,\
CD4_DNASE_IMPORT,\
GM06990_DNASE_IMPORT,\
Nessie_NG_STD_2_ctcf_ren_BR1

TARGET=\
CD4_H2AZ,\
CD4_H2BK5me1,\
CD4_H3K27me1,\
CD4_H3K27me2,\
CD4_H3K27me3,\
CD4_H3K36me1,\
CD4_H3K36me3,\
CD4_H3K4me1,\
CD4_H3K4me2,\
CD4_H3K4me3,\
CD4_H3K79me1,\
CD4_H3K79me2,\
CD4_H3K79me3,\
CD4_H3K9me1,\
CD4_H3K9me2,\
CD4_H3K9me3,\
CD4_H3R2me1,\
CD4_H3R2me2,\
CD4_H4K20me1,\
CD4_H4K20me3,\
CD4_H4R3me2,\
CD4_PolII,\
Wiggle_H3K27me3,\
Wiggle_H3K36me3,\
Wiggle_H3K4me3,\
Wiggle_H3K79me3,\
Wiggle_H3K9me3,\
Wiggle_H4K20me3

OUTDIR=/lustre/scratch1/ensembl/graef/RegBuild/v47

# first dump annotated features from database 
# for further processing

bsub -o "$OUTDIR/RegulatoryBuild_$FOCUS.log" -J RegBuild_Dump \
$EFG_SRC/scripts/build_regulatory_features.pl \
    -host $HOST \
    -port $PORT \
    -user $USER \
    -pass $PASS \
    -dbname homo_sapiens_funcgen_47_36i \
    -data_version 46_36h \
    -outdir $OUTDIR \
    -dump \
    -focus $FOCUS \
    -target $TARGET \
    $ARGS

# build regulatory features and dump to file when dumping 
# features has been sucessfully completed

bsub -o "$OUTDIR/RegulatoryBuild_$FOCUS.log" -J RegBuild_[1-25] \
    -w 'done(RegBuild_Dump)' \
    $EFG_SRC/scripts/build_regulatory_features.pl \
    -host $HOST \
    -port $PORT \
    -user $USER \
    -pass $PASS \
    -dbname homo_sapiens_funcgen_47_36i \
    -data_version 46_36h \
    -focus $FOCUS \
    -target $TARGET \
    -outdir $OUTDIR \
    -gene_signature \
    -dump_features \
    $ARGS

#    -seq_name 20 \
