#!/bin/sh

#Test for ENV_NAME here?
#Don't actually need to be in efg.env to run this
#Altho some default are inherited from the env.

usage='run_import_BED.sh $USER $PASS /your/BED/file/input_dir [ file_suffix ]\n\nInput dir will be used as the set suffix and should reflect the project name\n\te.g. /your/root/alignment/dir/homo_sapiens/GRCh37/ENCODE_Duke\nThe input files should be named CELLTYPE_FEATURETYPE.*.SUFFIX where the default SUFFIX=samse.bed.gz\n\te.g. GM12878_PolII.samse.bed.gz\nOptional file_suffix arg redefines input file SUFFIX'

USER=$1
shift

if [ ! "$USER" ]; then
	echo -e "ERROR:\tYou must provide a user argument.\n\n$usage"
	exit
fi


PASS=$1
shift

if [ ! "$PASS" ]; then
	echo -e "ERROR:\tYou must provide a password argument.\n\n$usage"
	exit
fi

input_dir=$1
shift

if [[ ! -d "$input_dir" ]]; then
	echo -e "ERROR:\tYou must provide a valid input_dir argument.\n\n$usage"
	exit
fi

suffix=$1
shift
suffix=${suffix:='samse.bed.gz'};


#Get full path
bed_dir=$(readlink -f $bed_dir)


#Vendor is only appropriate for input with vendor specific formats
#Use generic technology tag otherwise e.g. SEQUENCING
vendor='SEQUENCING'
species='homo_sapiens'
data_dir='/lustre/scratch103/ensembl/funcgen'
dbname="dev_${species}_funcgen_59_37d"


#Validate here by asking experiment name?
#This is not really experiment name
#more experiment_group/set InputSet suffix


project=$(echo $input_dir  | sed 's/\/$//')
project=$(echo $project | sed 's/.*\///')

for filepath in $(ls $input_dir/*$suffix); do

	file=$(echo $filepath | sed 's/.*\///')
		
	input_set=$(echo $file | sed 's/\.samse.*//')
	cell_type=$(echo $input_set | sed 's/_.*//')
	feature_type=$(echo $input_set | sed -r "s/${cell_type}_//")
	feature_type=$(echo $feature_type | sed 's/_.*//')
	input_set="${input_set}_${project}"

	
	if [[ $feature_type != WCE ]]; then
		echo -e "Project:\t$project"
		echo -e "InputSet:\t$input_set"
		echo -e "CellType:\t$cell_type"
		echo -e "FeatureType:\t$feature_type"
		
#Now we are preparing files first, we need to submit this to the farm
#otherwise this loop my take a long time and a lot of CPU on the head nodes

		job_name="parse_and_import_${input_set}"
		output_dir="${data_dir}/output/${dbname}/${vendor}/${input_set}"
        #Create here as Importer logs try to access this prior to creation

		
		if [[ ! -d "$output_dir" ]]; then
			mkdir -p $output_dir
		fi

#Could do with -R " rusage [tmp=filesize] "
#This would stop multiple job running on the same node and failing due to lack of tmp space




		bsub_cmd="bsub -q long -J $job_name -o ${output_dir}/${job_name}.out -e ${output_dir}/${job_name}.err"
#This will still log to the normal file, only LSF out and err output will be caught here

#Will this still cause /tmp space problems if the are run on the same host?
#These are quite often all submitted to the same host if possible max 4jobs/host
#get size of file + 1GB
#-R"select[tmp>20G] rusage[tmp=20G]" ?


		$bsub_cmd time perl $EFG_SRC/scripts/import/parse_and_import.pl\
			-name  $input_set\
			-format SEQUENCING\
			-parser Bed\
			-vendor $vendor\
			-location Hinxton\
			-contact njohnson@ebi.ac.uk\
			-group efg\
			-species $species\
			-input_set $input_set\
			-input_feature_class result\
			-registry_host ens-livemirror\
			-registry_port 3306\
			-registry_user ensro\
			-assembly 37\
			-host ens-genomics1\
			-port 3306\
			-dbname $dbname\
			-user $USER\
			-pass $PASS\
			-cell_type $cell_type\
			-feature_type $feature_type\
			-feature_analysis bwa_samse\
			-recover\
			-result_files $input_dir/$file\
			-data_root $data_dir\
			-farm\
             $@ 
#"

#-slice 1\
#This used to work with test sub slices but now complains when trying to set the slices???



#removed -tee to prevent redundancy between log and out file
#	-skip_slices 1
#exit
	fi
done

#Not we never project reads between assemblies here
#Run the mapping pipeline first!

	#-result_files $*      



#Really need to set IMPORTED status on input_subset for this file
#Can't do this currently as we run batched farm jobs
