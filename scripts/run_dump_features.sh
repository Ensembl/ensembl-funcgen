#!/bin/sh

USER=$1
shift

PASS=$1
shift

#test for user!

if [ ! $USER ]; then
	echo "You must specific a user argument"
	exit
fi

if [ $PASS ]; then
	PASS="-pass $PASS"
fi


#set_names="K562_SP1_ENCODE_Hudsonalpha_SWEMBL_R015 HepG2_USF1_ENCODE_Hudsonalpha_SWEMBL_R015"
#set_names="RegulatoryFeatures:MultiCell" # RegulatoryFeatures:K562 RegulatoryFeatures:CD4 RegulatoryFeatures:IMR90  RegulatoryFeatures:GM06990 RegulatoryFeatures:GM12878 RegulatoryFeatures:H1ESC RegulatoryFeatures:H1ESC RegulatoryFeatures:HepG2 RegulatoryFeatures:NHEK RegulatoryFeatures:HUVEC"

dbport=3306
dbname='homo_sapiens_funcgen_62_37g'
dbhost='your_db_host'
format='bed'
feature_type='ProbeFeature'

#set_names="RegulatoryFeatures:MultiCell RegulatoryFeatures:ES RegulatoryFeatures:ESHyb RegulatoryFeatures:MEF RegulatoryFeatures:NPC"
#dbname='mus_musculus_funcgen_62_37o'

out_dir="your_out_dir"

if [ ! -d $out_dir ]; then
	mkdir -m 775 -p $out_dir
fi


#for set_name in $set_names; do
#make outdir for lsf out
#

#dump_name=$set_names
dump_name='HG-U133A'
#dump_name='sets'
#dump_params="-feature_sets $set_names -dump_name $dump_name"
dump_params="-array $dump_name -vendor AFFY"

#cdb host is hard coded to staging, needs changing
	bsub_cmd="bsub -q long -J dump_gff_${set_name} -o ${out_dir}/dump_gff_${dump_name}.out -e ${out_dir}/dump_gff_${dump_name}.err -R\"select[mem>6000] rusage[mem=6000]\" -M 6000000"
	#MultiCell tends to run out of memory!

	job_cmd="$EFG_SRC/scripts/export/dump_features.pl\
	-port   $dbport\
    -format $format\
	-dbhost $dbhost\
	-dbname $dbname\
	$dump_params\
	-outdir $out_dir \
	-user $USER\
	$PASS\
	$@	"


#-slices 1

	echo "$bsub_cmd $job_cmd"


	#Have to echo through bash to stop weird splitting of -R option
	#echo "$bsub_cmd  $job_cmd" | bash
	$job_cmd
#done



exit;
