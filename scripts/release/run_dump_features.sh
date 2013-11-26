#!/bin/sh

# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


USER=$1
shift

if [ ! $USER ]; then
	echo "You must specific a user argument"
	exit
fi


format='GFF'
release=72

#only required if we want to dump specific sets i.e. not the complete reg build set
#fset_names="RegulatoryFeatures:MultiCell RegulatoryFeatures:ES RegulatoryFeatures:ESHyb RegulatoryFeatures:MEF RegulatoryFeatures:NPC RegulatoryFeatures:MEL"

#dbname="homo_sapiens_funcgen_${release}_37"
#dbhost='ens-staging1'

dbname="mus_musculus_funcgen_${release}_38"
dbhost='ens-staging2'



#feature_type='AnnotatedFeature'
#feature_type='RegulatoryFeature'
feature_type='MotifFeature'
#feature_type='SegmentationFeature'

fset_names="${feature_type}s" #Special set name to dump all Annotated/MotifFeatures in build
#Does not yet support RegulatoryFeatures as these seem to get dumped to a merged file with no cell type info
#and merge does not handle this yet anyway



# MySQL query to get the correct $fset_names:
# select replace(group_concat(name), ",", " ") from feature_set where type = 'regulatory' ;

# Mouse:
#fset_names="RegulatoryFeatures:ESHyb RegulatoryFeatures:ES RegulatoryFeatures:NPC RegulatoryFeatures:MultiCell RegulatoryFeatures:MEF RegulatoryFeatures:MEL"

# Human:
# Regulatory
#fset_names="RegulatoryFeatures:NHEK RegulatoryFeatures:K562 RegulatoryFeatures:GM06990 RegulatoryFeatures:MultiCell RegulatoryFeatures:IMR90 RegulatoryFeatures:HSMM RegulatoryFeatures:HepG2 RegulatoryFeatures:NH-A RegulatoryFeatures:HeLa-S3 RegulatoryFeatures:CD4 RegulatoryFeatures:HUVEC RegulatoryFeatures:HMEC RegulatoryFeatures:H1ESC RegulatoryFeatures:GM12878";
# Segmentation
#fset_names="Segmentation:HUVEC Segmentation:K562 Segmentation:GM12878 Segmentation:H1ESC Segmentation:HeLa-S3 Segmentation:HepG2";
#Should do this for all, so we don't have to specify -feature_sets

dnadbhost=$dbhost
dnadbuser=$USER

root_dir='/lustre/scratch109/ensembl/funcgen/output'
out_dir="${root_dir}/${dbname}/dumps/${format}/${feature_type}"


#Really need to hoik out the job definition/submission code to a Hive module, which populates the input_ids accordingly
#Do we really need a Hive for this?

#These also need to have some DB load monitoring applied to the resource spec

if [[ $USER != merge ]]; then
    


	#Only do this loop for reg feat sets
	#i.e. we want separate dumps
	#remove for loop for AnnotatedFeature dumps i.e. merged dump

	#Could get script to automatically define rsets based on feature class but
	#Loop here to avoid having to do it in the script for local jobs


    for set_names in $fset_names; do

        dump_name="-dump_name $set_names"
        dump_params="-feature_sets $set_names $dump_name"

        #dump_params="-result_sets $set_names -window_size $wsize"
        #dump_params="-array $dump_name -vendor AFFY"
		
        bin_dir='/nfs/users/nfs_t/tj1/'
			
        if [ ! -d $out_dir ]; then
            mkdir -m 775 -p $out_dir
        fi

        job_cmd="$SRC/ensembl-functgenomics/scripts/export/dump_features.pl\
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
	
else # USER = merge
	
	#MERGE dumps

	species=$(echo $dbname | sed 's/_funcgen.*//')
	ftp_dir="${root_dir}/ftp/release-${release}/regulation/$species"

	#start_dir=$PWD

	for fset in $fset_names; do
      file_prefix=$(echo $fset | sed s'/\:/_/')
      dump_dir="${out_dir}/${fset}"
      echo -e "\nMerging $file_prefix $feature_type dumps\ncd'ing to $dump_dir"
      cd $dump_dir
		
	    #remove old merged dump first
		#so we don't unzip it below
		if [[ -e  ${file_prefix}.gff.gz ]]; then
			rm -f ${file_prefix}.gff.gz
		fi


		#Unzip files if required		
		files=$(ls ${file_prefix}.*.gff.gz 2>/dev/null) 

		if [[ $files ]]; then
			gunzip $files
		fi
		
		
		cmd="cat ${file_prefix}.*.gff | gzip > ${file_prefix}.gff.gz"
		echo $cmd
		eval $cmd
		#can't simply call $cmd here for some reason as cat's to stdout
		#hence need eval
		


		#Never remove source dumps in case cat fails

		if [[ ! -d $ftp_dir ]]; then
        mkdir -p $ftp_dir
		fi
    
		cmd="ln -sf $PWD/${file_prefix}.gff.gz ${ftp_dir}/${file_prefix}.gff.gz"
		echo $cmd
		$cmd

				
    done


	#Now generate the md5s in the ftp dir
	#Put this in a separate script, generate_md5s_by_file_suffix.sh
	#Can then re-use for data_files

	cd $ftp_dir
	checksum_file="${format}_hex_checksums.MD5"

	if [ -e $checksum_file ]; then
		rm -f $checksum_file
	fi

	for f in $(ls ./); do

		#This

		if [ $f != README ]; then

			#when f=RegulatoryFeatures_NH-A.gff.gz

		    #Warning: Use of "-A" without parentheses is ambiguous at -e line 1.
			#syntax error at -e line 1, near "RegulatoryFeatures_NH-A"
			#Execution of -e aborted due to compilation errors.

			#Also .'s in file name are being missed in print?

	  #do we need an eval here?

			perl -MDigest::MD5 -e "open(FILE, './$f') or die(\"Cannot open file for digest:\t$f\"); print Digest::MD5->new->addfile(*FILE)->hexdigest.\" $f\n\";" >> $checksum_file
		fi

	done

	echo "DON'T FORGET TO COPY THE README FROM THE PREVIOUS RELEASE!"
	echo "Are the permissions set correctly? 775?"
fi
