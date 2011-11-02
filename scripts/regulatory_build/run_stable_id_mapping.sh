#!/bin/sh


USER=$1
shift
PASS=$1
shift

if [ ! $USER ] || [ ! $PASS ]; then
 echo 'Usage: run_stable_id_mapping.sh $USER $PASS [ -other args ]'; 
 exit 1;
fi

species='homo_sapiens'
#species='mus_musculus'
data_dir='/lustre/scratch103/ensembl/funcgen/output/'
ndbname="devnew_${species}_funcgen_65_37"
out_dir="${data_dir}/${ndbname}/stable_id_mapping"
ndbhost='ens-genomics2'
odbname=$ndbname
odbhost=$ndbhost

if [ ! -d $out_dir ]; then
	mkdir -p $out_dir
fi

bsub_cmd="bsub -J ${ndbname}_stable_id_mapping -q yesterday -o $out_dir/stable_id_mapper.%J.$$.out -e $out_dir/stable_id_mapper.%J.$$.err"

job_cmd="perl -w $EFG_SRC/scripts/regulatory_build/stable_id_mapper.pl\
 	-ndbname $ndbname \
 	-nhost   $ndbhost \
 	-nuser   $USER \
   	-npass   $PASS \
    -ouser   ensro\
    -ohost   $odbhost
    -odbname $odbname\
    -old_fset_name RegulatoryFeatures:MultiCell_v11\
	-new_fset_name RegulatoryFeatures:MultiCell\
	-out_dir  $out_dir\
    -species  $species\
  -log_file ${out_dir}/stable_id_mapper.$$.log
 $@
"


echo "$job_cmd"

$bsub_cmd $job_cmd



exit;

#Currently have to -clobber to overcome human Y par issues 

#	-old_assembly NCBI36\
#	-new_assembly GRCh37\
#	$@

#-next_stable_id
#-slices
#-assign_all_nulls

#stable ID mapping bug may be due to build using non-ref seq_regions which have features assembly mapped from the ref seq?
#The new seq_regs do no have a comparable old seq_reg to map from, hence they are simply assigned a new stable ID



#This is all now done in the stable_id_mapper.pl

#Test we have no NULL stable IDs for MultiCell

#select fs.name, sr.seq_region_id, sr.name, count(distinct rf.regulatory_feature_id) from seq_region sr, regulatory_feature rf, feature_set fs where fs.feature_set_id=rf.feature_set_id and stable_id is NULL and fs.name='RegulatoryFeatures:MultiCell' and rf.seq_region_id=sr.seq_region_id group by fs.feature_set_id, seq_region_id;


#HC sql here
#select count(*) from feature_set fs, feature_set fs1, regulatory_feature rf, regulatory_feature rf1 where rf1.seq_region_id=rf.seq_region_id and rf1.seq_region_start=rf.seq_region_start and rf1.seq_region_end=rf.seq_region_id and rf.feature_set_id=fs.feature_set_id and fs.name rlike 'RegulatoryFeatures:MultiCell_v[0-9]+' and rf1.feature_set_id=fs1.feature_set_id and fs1.name='RegulatoryFeatures:MultiCell' and rf.stable_id!=rf1.stable_id

#update regulatory_feature rf, regulatory_feature rf1 left join feature_set fs1 on rf1.feature_set_id=fs1.feature_set_id set rf1.stable_id=rf.stable_id where rf.seq_region_start=rf1.seq_region_start and rf.seq_region_id=rf1.seq_region_id and rf.feature_set_id=(select feature_set_id from feature_set where name='RegulatoryFeatures:MultiCell') and fs1.name not rlike ".*_v[0-9]+" and fs1.name!='RegulatoryFeatures:MultiCell';


#Test we have no NULL stable IDs for any set
#select fs.name, count(distinct rf.regulatory_feature_id) from regulatory_feature rf, feature_set fs where fs.feature_set_id=rf.feature_set_id and stable_id is NULL group by fs.feature_set_id;


