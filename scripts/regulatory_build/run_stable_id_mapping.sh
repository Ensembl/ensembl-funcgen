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



#Currently have to -clobber to overcome human Y par issues 

#	-old_assembly NCBI36\
#	-new_assembly GRCh37\
#	$@

#-next_stable_id
#-slices
#-assign_all_nulls

#Need to add this to the stable_id_mapper
update regulatory_feature rf, regulatory_feature rf1 left join feature_set fs1 on rf1.feature_set_id=fs1.feature_set_id set rf1.stable_id=rf.stable_id where rf.seq_region_start=rf1.seq_region_start and rf.seq_region_id=rf1.seq_region_id and rf.feature_set_id=(select feature_set_id from feature_set where name='RegulatoryFeatures:MultiCell') and fs1.name not rlike ".*_v[0-9]+" and fs1.name!='RegulatoryFeatures:MultiCell';


#Test we have no NULL stable IDs for any set
select fs.name, count(distinct rf.regulatory_feature_id) from regulatory_feature rf, feature_set fs where fs.feature_set_id=rf.feature_set_id and stable_id is NULL group by fs.feature_set_id;


