#!/bin/sh

if [ $# -lt 4 ]; then

	echo "Usage: $0 <host> <port> <user> <password> <mode>"
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
MODE=$1
shift
ARGS=$*

DATA_VERSION='51_36m'
DBNAME='homo_sapiens_funcgen_'${DATA_VERSION}

OUTDIR='/lustre/work1/ensembl/graef/RegBuild/v51'

FOCUS=\
CD4_CTCF_e51,\
CD4_DNASE_IMPORT,\
GM06990_DNASE_IMPORT,\
Nessie_NG_STD_2_ctcf_ren_BR1

TARGET=\
CD4_H2AZ_e51,\
CD4_H2BK5me1_e51,\
CD4_H3K27me1_e51,\
CD4_H3K27me2_e51,\
CD4_H3K27me3_e51,\
CD4_H3K36me1_e51,\
CD4_H3K36me3_e51,\
CD4_H3K4me1_e51,\
CD4_H3K4me2_e51,\
CD4_H3K4me3_e51,\
CD4_H3K79me1_e51,\
CD4_H3K79me2_e51,\
CD4_H3K79me3_e51,\
CD4_H3K9me1_e51,\
CD4_H3K9me2_e51,\
CD4_H3K9me3_e51,\
CD4_H3R2me1_e51,\
CD4_H3R2me2_e51,\
CD4_H4K20me1_e51,\
CD4_H4K20me3_e51,\
CD4_H4R3me2_e51,\
CD4_PolII_e51,\
Wiggle_H3K27me3,\
Wiggle_H3K36me3,\
Wiggle_H3K4me3,\
Wiggle_H3K79me3,\
Wiggle_H3K9me3,\
Wiggle_H4K20me3,\
CD4_H2AK5ac_e51,\
CD4_H2AK9ac_e51,\
CD4_H2BK120ac_e51,\
CD4_H2BK12ac_e51,\
CD4_H2BK20ac_e51,\
CD4_H2BK5ac_e51,\
CD4_H3K14ac_e51,\
CD4_H3K18ac_e51,\
CD4_H3K23ac_e51,\
CD4_H3K27ac_e51,\
CD4_H3K36ac_e51,\
CD4_H3K4ac_e51,\
CD4_H3K9ac_e51,\
CD4_H4K12ac_e51,\
CD4_H4K16ac_e51,\
CD4_H4K5ac_e51,\
CD4_H4K8ac_e51,\
CD4_H4K91ac_e51

case "$MODE" in 
    "dump")


        # first dump annotated features from database 
        # for further processing

        bsub -o "$OUTDIR/RegulatoryBuild_afDump_$FOCUS.log" -J RegBuild_afDump \
        $EFG_SRC/scripts/regulatory_build/build_regulatory_features.pl \
            -host $HOST \
            -port $PORT \
            -user $USER \
            -pass $PASS \
            -dbname $DBNAME \
            -data_version $DATA_VERSION \
            -outdir $OUTDIR \
            -dump \
            -focus $FOCUS \
            -target $TARGET \
            $ARGS
        ;;
    
    * )

        # build regulatory features and dump to file when dumping 
        # features has been sucessfully completed

        #

        bsub -o "$OUTDIR/RegulatoryBuild_$FOCUS.log" -J RegBuild_[1-25] \
            -w 'done("RegBuild_afDump")' \
        $EFG_SRC/scripts/regulatory_build/build_regulatory_features.pl \
            -host $HOST \
            -port $PORT \
            -user $USER \
            -pass $PASS \
            -dbname $DBNAME \
            -data_version $DATA_VERSION \
            -focus $FOCUS \
            -target $TARGET \
            -outdir $OUTDIR \
            -dump_features \
            -clobber \
            -stats \
            $ARGS
        
           # -gene_signature \
           # -write_features \
           # -seq_name 10 \
       ;;
esac
          
