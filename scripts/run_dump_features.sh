#!/bin/sh

USER=$1
shift

if [ ! $USER ]; then
	echo "You must specific a user argument"
	exit
fi



#set_names="K562_SP1_ENCODE_Hudsonalpha_SWEMBL_R015 HepG2_USF1_ENCODE_Hudsonalpha_SWEMBL_R015"
set_names=""



#rset_names="RegulatoryFeatures:MultiCell RegulatoryFeatures:K562 RegulatoryFeatures:CD4 RegulatoryFeatures:IMR90  RegulatoryFeatures:GM06990 RegulatoryFeatures:GM12878 RegulatoryFeatures:H1ESC RegulatoryFeatures:H1ESC RegulatoryFeatures:HepG2 RegulatoryFeatures:NHEK RegulatoryFeatures:HUVEC"


#dbname='homo_sapiens_funcgen_63_37'
#dbhost='ens-staging'
#format='bigwig'
format='GFF'
#feature_type='ProbeFeature'
#feature_type='ResultFeature'
feature_type='RegulatoryFeature'
#feature_type='AnnotatedFeature'



rset_names="RegulatoryFeatures:MultiCell RegulatoryFeatures:ES RegulatoryFeatures:ESHyb RegulatoryFeatures:MEF RegulatoryFeatures:NPC"
#rset_names='AnnotatedFeatures' #Special set name to dump all supportingt sets
dbname='mus_musculus_funcgen_63_37'
dbhost='ens-staging2'

dnadbhost=$dbhost
dnadbuser=$USER


#change output to RegulatoryFeature.ctype rather than RegulatoryFeatures_ctype?

#set_names="IMR90_DNase1_EpiRoadmap_Uw_GSE18927"
#set_names="AnnotatedFeatures"

if [[ $USER != merge ]];then


	#Only do this loop for reg feat sets
	#i.e. we want separate dumps
	#remove for loop for AnnotatedFeature dumps i.e. merged dump

	for set_names in $rset_names; do



dump_name="-dump_name $set_names"
#dump_name='-dump_name HG-U133A'

dump_params="-feature_sets $set_names $dump_name"
# [30, 65, 130, 260, 450, 648, 950, 1296];
#wsizes="30 65 130 260 450 648 950" # 1296"
#wsizes="1296"
#for wsize in $wsizes; do


#dump_params="-result_sets $set_names -window_size $wsize"
#dump_params="-array $dump_name -vendor AFFY"

bin_dir='/nfs/users/nfs_n/nj1/'
out_dir="/lustre/scratch103/ensembl/funcgen/output/${dbname}/dumps/${format}/${feature_type}"

if [ ! -d $out_dir ]; then
	mkdir -m 775 -p $out_dir
fi

	job_cmd="$EFG_SRC/scripts/export/dump_features.pl\
	-port 3306\
    -format $format\
	-dbhost $dbhost\
	-dbname $dbname\
    -dnadb_host $dnadbhost\
    -dnadb_user $dnadbuser\
	$dump_params\
	-out_root $out_dir \
    -bin_dir $bin_dir\
	-user $USER\
	$@	"


#Add these on cmd line to avoid including them by mistake
#-farm
#-post_process
#-force_local
#-slices chromosome:GRCh37:X:1:421628\

	$job_cmd
done


else

	#MERGE RegulatoryFeature dumps

	if [[ $feature_type != RegulatoryFeature ]]; then
		echo "Can only merge RegulatoryFeature sets at present"
	fi

species=$(echo $dbname | sed 's/_funcgen.*//')
ftp_dir=/lustre/scratch103/ensembl/funcgen/output/ftp/current_functional_genomics/$species


start_dir=$PWD
cd $out_dir

for rset in $rset_names; do
	file_prefix=$(echo $rset | sed s'/\:/_/')

	gunzip ${file_prefix}*gff.gz
	rm -f ${file_prefix}.gff

	cat ${file_prefix}*gff > ${file_prefix}.gff

	if [[ ! -d $ftp_dir ]]; then
		mkdir -p $ftp_dir
	fi

	#Add MD5s here and test?

	ln -s $PWD/${file_prefix}.gff ${ftp_dir}/${file_prefix}.gff
done

echo "DON'T FORGET THE READMEs!"

cd $start_dir

fi


#Are the permissions set correctly? 775?
