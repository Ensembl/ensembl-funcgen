#!/bin/sh



HOST=your_host
PORT=3306
USER=your_user
PASS=$1
shift
ARGS=$*

if [ ! $PASS ]; then
	echo "You must supply a password argument"
fi


DBNAME=your_db_name
OUTDIR="$EFG_DATA/reg_build/${DBNAME}"

if [ ! -d $OUTDIR ]; then
	mkdir -p $OUTDIR
fi


# Old params superceded by -use_tracking
#Comma or space separated list of focus and attribute sets
#Sets with spaces in names must be comma separated and quoted appropriate
#FOCUS="\
#FOCUS_SET1,\
#FOCUS_SET2,\
#FOCUS_SET2"

#ATTRIBUTE sets must only reflect CellTypes present in FOCUS sets
#ATTRIBUTE="\
#ATTR_SET1,\
#ATTR_SET2,\
#ATTR_SET3,\
#ATTR_SET4,\
#ATTR_SET5,\
#ATTR_SET6,\
#ATTR_SET7,\
#ATTR_SET8"

VERSION=your_incremented_version_number

$EFG_SRC/scripts/regulatory_build/build_regulatory_features.pl\
    -host $HOST \
	-port $PORT \
	-user $USER \
	-pass $PASS \
	-dbname $DBNAME \
	-use_tracking\
	-outdir $OUTDIR \
	-version $VERSION\
	-write_features\
	-bsub_mem 10000\
	-rollback\
	-tee\
	-logfile $OUTDIR/RegulatoryBuild.log\
	$ARGS

#Only required if current assembly is not on ensembldb
#or if -gene_signature is specified
#  -dnadb_name your_dnadb_name\
#  -dnadb_host your_dnadb_host\

# Old params superceded by -use_tracking
#	-focus_sets $FOCUS \
#	-attribute_sets $ATTRIBUTE \


#For more help use:
#build_regulatory_features.pl -help
#or
#perldoc build_regulatory_features.pl
