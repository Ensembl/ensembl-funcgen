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
# LICENCE:
# This code is distributed under an Apache style licence. Please see
# http://www.ensembl.org/info/about/code_licence.html for details.
#
# CONTACT
# Please post comments/questions to the Ensembl development list
# <ensembl-dev@ebi.ac.uk>



echo ':::: Welcome to the eFG array mapping environment'

#export trap_exit=1
#Also change so we don't need absolute path as we can source from the same directory 

. $EFG_SRC/scripts/environments/pipeline.env


export ENV_NAME='arrays'
#env colour is green
export PS1_COLOUR=32
ARRAYS_CONFIG=${ARRAYS_CONFIG:=$EFG_SRC/scripts/environments/arrays.config}
echo ":: Sourcing ARRAYS_CONFIG: $ARRAYS_CONFIG"
. $ARRAYS_CONFIG


### DONE
# DONE Error catch properly > Execute
# DONE implement 'Execute' for critical commands, see func.sh and affy.env Execute
# DONE Add full setup for each rule in config. Each array in ImportArrays, Each format in Probe(Transcript)Align
# DONE DumpSoftMaskedSeq should handle striping on GENOMICSEQS dir.
# DONE Collapses each array format separately merging to one NR_FASTA. Different ARRAY_DESIGN should be run separately.
# DONE enable multiple arrays.list files for each array format
# DONE implement ARRAY_ENV to reduce the amount of vars set e.g. do not extend PERL5LIB ?
# DONE Enable parallelised ImportArrays using dynamic bsub config
# DONE Base ImportArrays bsub -R config on size of files in each format group
# DONE Move All set up to init sub, clean all vars when sourcing, then call init after species_VERSION.arrays has finished set up.  This will stop vars being inherited from previous envs.
# DONE Change ARRAY_DESIGN check to use VALID_DESIGNS array, can we grep this?
# DONE Ask if we want to dump GENOMICSEQS/TRANSCRIPTSEQS if not set. Also test files
# DONE Change BackUpFile to use back up dir if set.
# DONE Change OUTPUT_DIR to WORK_DIR and genera
# DONE Source .pipeline
# DONE Use awk to validate FASTA
# DONE Sort BuildFasta so we don't mis any sub string arrays when greping?
# DONE No way of doing this ? As list always contains full file name which could zip name or fasta name if it wasn't zipped
# DONE To get around this, we need to unzip all fastas and add unzipped name to list, then match on exact name
# DONE This depends on unzipped file being the same as the zipped file s/zip|gz//
# DONE RollbackArrays now handles all previous 'Clean'/rollback functionality via rollback_array.pl/Helper.pm
# DONE Remove DEFAULT_ARRAY_FORMATS, now uses ARRAYS_HOME or ARRAY_FORMATS in dbname.arrays
# DONE MakeHugememWorkdir
# DONE Remove SPECIES_COMMON? Now tr to uc SPECIES
# DONE BuildFastas shoudl always build from scratch to avoid reusing badly formated data.
# DONE Implement QUEUE_MANAGER to set "local" modes, also use env var to set BatchQueue.pm QUEUE_MANAGER Need to default to LSF if not set
# DONE Auto generate BatchQueue.pm config
# DONE Handle QUEUE_MANAGER in submitJob/checkJob, implemented in BatchQueue.pm
# DONE Remove SPECIES_COMMON and just uc SPECIES to get the arrays_dir
# DONE Move all vars to separate file and just have defaults in here.
# DONE Distribute/ArchiveData
# DONE Need to account for align_types in SetUpPipeline. Implement ALIGN_TYPES similar to ARRAY_FORMATS

### TO DO LIST ###
# Add exit code handling, 101 for miscleaneous.
# add ARRAY_FORMAT to arrays.list file and BuildFasta method, so custom list works with multiple formats
# Catch ImportArrays $? error and also archive old output files, so whe can test for the rpesenc of new output files when cat'ing
# Add ISKIPLIST/REGEX/FIELD to config to enable skipping of certain non-experimental probes
# Move all interactive stuff to Init to enable full automation of RunAlignments
# Set vars for all interactive questions to allow full non-interactive automation, need to Execute here to error catch
# Make batch size dynamic
# Create CreateAlignIDs RunnableDB? And runnable DBs for all array wide clean steps?
# Check if we are handling haplotype regions in dump, did I fix the memory problem?
# Use the db fasta cache?
# Test DB connections in Init
# Grep for mysql reserved characters in DBNAME in _Init? Or just capture mysql failure better.
# Use default seq file based one validated build?
# Need to reset jobs align jobs when we CleanArrays?
# Roll back by chunk(This would require using/recreating exonerates chunk function). 
# Currently no way of rolling back output for individual chunk
#   So when we re-run we may have duplicate Unmapped objects and features depending on where the failure happened
#   This can result in anomalous Transcript annotations due to duplicate features which are not caught by probe2transcript
# Create CleanUpProbeAligns: We are currently storing unmapped objects for *_ProbeTranscriptAlign which may be spurious
# due to valid alignments in *_ProbeAlign. And vice versa, we may have exceeded hits in ProbeAlign, but ProbeTranscriptAlign
# will throw these away, so will not be aware that it is a promiscuous probe. Will have to update Helper::rollback_ArrayChip. We could integrate the Genomic and Transcript alignment into the same job, running one after the other?
# This would enable us to use the same caches and therefore would remove clean up step.
# GenerateCDF
# Man pages for functions?
# Make sure we have all the functionality friom these scripts $SRC/ensembl-personal/npj/oligo_mapping/scripts
# array_FORMAT.names can empty if ImprtArrays fell over, we should warn about this in Rollback. Can also become empty after ImportArrays somehow?
# Add LSF DB throttling in SubmitAlign?
# Set up individual AlignWaits for each analysis, then we can kick off RunTranscriptXrefs straight away?
# This would require threading or integration of probe2transcript into a RunnableDB.
# Same here for CreateAlignIDs step and Import to Align transition
# Make probe_per_chunk and batch size dependant on size of TRANSCRIPT/GENOMICSEQS file and array_nr for each format
# This is also how we can set the -R separately for each of the CONFIG by setting \$${FORMAT}_CONFIG_VAR before submitting
# and having $"FORMAT"_CONFIG_VAR in config for each FORMAT 
# THIS IS WHY WE HAVE TO DO IT IN STEPS AS EVEN IF WE COULD INCLUDE SETTING THE ENV VARS IN IMPORTARRAYS
# THESE WOULD NOT BE PROPOGATED TO THE ENV OF THE PIPELINE(just the import job env)
# Same for probe2transcript submission? Altho these could be calculated in CreateAlignIDs step.
# Implement NetAffx script to check for new arrays?
# Allow adding format support to probe2transcript without editing code directly
# PROBE_ALIGN_RESOURCE: Some AFFY_UTR_ProbeAlign jobs are running out of memory for chimp.
# genomics seqs file should be dumped with assembly version, not just schema version!
# Add func to compare config file with config_template to identify whether there are any missed updates
# Load rules/analyses in line to avoid warnings during set up?

# ISSUE: we cannot re-do transcript mapping without reconstituting the nr fasta file with dbIDs
# Cannot assume we will always have the original file, so we need to alter ImportArrays to run with the original sources files, but simply use the DB to retrieve dbIDs, rather than storing the probes.  We could extend this to allow incremental updates of new arrays rather than having to clear the DB and re import everything.
# ISSUE: Some alt-trans may give rise to duplicate ProbeTranscriptAlign mappings.  This should not currently effect probve2transcript as we are simply count the distinct probes.  But we should really count a key ENSG|:probe_id:genomic_start:genomic_stop.  Then we can count multiple mappings of the same probe


#DOCS
#Need to define format and array specific config in Analysis/Config/ImportArrays
#ARRAY_PARAMS needs to match the content of the file, not the name of the file as this can differ.
#Need to set up format dirs in species home reflecting these array formats
# add doc to the top of this env about dir setup, and required code bases e.g. analysis, pipeline, core, efg
#Add docs to the VALID vars, about what you would need to do to add an an extra valid var, and a summary description of how they are used


### Some convinient aliases
# Some also set in pipeline.env
alias arraysdir='cd $ARRAYS_HOME' 




################################################################################
# Func      : _InitEnv(
# Desc      : Sets up all config defined variables. Should only be called from
#             instance config file (species_VERSION.arrays) e.g. mouse_51.arrays 
# Args [n]  : 
# Return    : none
# Exception : 
################################################################################

_InitEnv(){

	#Set all generic pipeline stuff e.g. DB vars, PATH/PERL5LIB etc
	_InitPipelineEnv

	
	#Do a getopts here?
	#How do we make this pass to _InitPipelineEnv?
	#use \?) to build other params string?
	#Will this handle args properly? May have to write BuildUnkownOptions function?
	#This will grab params and reset $OPTIND accordingly
	#This is only for dump/use seqs with matching schema_builds
	#create db



	echo ":: Setting config for $SPECIES array mapping"

	if [ $warn ]
	then
		echo 'NOTE: NIMBLEGEN_TILING will never go through ImportArrays or ProbeTranscriptAlign steps?'
	fi
	
	#ValidateBooleans can we write a method for this
	#Move this to pipeline.env?

	if [ $MULTI_SPECIES ]; then

		if [ $MULTI_SPECIES -ne 1 ]; then
		echo 'MULTI_SPECIES is a boolean variable, please specify 1 or omit for default of 0'
		else
		#change to script paramter
			export MULTI_SPECIES=' -multi_species '
		fi
	fi
	

	#Let's trim this and move the host and working dir to the window title if we can
	#Can we set these so they are dynamic dependent on the value of $DATA_HOME? 
	#i.e. to facilitate use of $HOME for hugemem I/O
	#Contains only array format data dirs or the GENOMIC or TRANSCRIPTSEQS dirs
	#These should also match config in Analysis/Config/ImportArrays
	#uc SPECIES_COMMON here?


	#TO DO Need to Check config vars here now we have moved them to a separate configfile.
	#CheckVariables
	


	export ARRAYS_HOME=${DATA_HOME}/$(echo $SPECIES | tr [a-z] [A-Z])

	export WORK_DIR=${DATA_HOME}/${DB_NAME}
	#This is used in pipeline alias workdir
	export BACKUP_DIR=${WORK_DIR}/backup
	export PIPELINE_OUT=${WORK_DIR}/align_out

	#Check some dirs are present
	error=$(CheckDirs $ARRAYS_HOME)

	#This can be put in pipeline.env, even thos it is dependant on $DATA_HOME which is only defined in the 
	#top level environment e.g. arrays or peaks

	if [ $? -ne 0 ]; then
		echo $error;
		array_home_subdir=$(echo $ARRAYS_HOME | sed -r "s'^${DATA_DIR}''")
		echo -e "Attempting to get data from archive:\t$ARCHIVE_DIR/$array_home_subdir"
				#Need to get the diff between DATA_DIR and DATA_HOME
		#Also need to define these properly?
		DistributeData $ARCHIVE_DIR/${array_home_subdir}/


		if [ $? -ne 0 ]; then
			echo "Failed to find input data, please create and populate \$ARRAYS_HOME. See ensembl-functgenomics/docs/array_mapping.txt for details"
		fi

	fi

	#We also want to check whether we have a version in $ARCHIVE_DIR and DistributeData if we have
	#We actually want to create this, then exit suggesting format subdirs need to be created


	#Create db output dir if not present
	MakeDirs $WORK_DIR

	## Test and dump seq files ##
	## Also for TRANSCRIPT ALIGN_TYPE do extra DNADB checks
	## Import transcripts as coord_system if not present

	#Actually version is not right test for transcript set
	#Validate version
	#version_db=$(echo $DB_NAME | grep /[a-z]+_funcgen_${VERSION
	#check standard location and match version

	#Can we add /data/blastdb paths to here, no standard naming :( !!!

	if [ $warn ]
	then
		echo "NOTE:: Separate these into a separate sub to avoid having lengthy seq_dumping happening in the init func?"
		echo "NOTE:: Or do we want this interactive stuff forced at the start, so we don't nohup this and then have it hanging?"
	fi


	#Set/validate formats/types
	SetArrayFormats $ARRAY_FORMATS
	SetAlignTypes $ALIGN_TYPES
	jids=

	for align_type in $ALIGN_TYPES
	do
	
		#Check/SetUp TRANSCRIPT dependant vars


		#Can remove this as this is generic for both seq dumps

		if [ $align_type = TRANSCRIPT ]
		then

			if [ $warn ]
			then
				echo 'NOTE: Also need to check DNADB vars if ALIGN_TYPE is TRANSCRIPT'
			fi

			CheckVariables DNADB_NAME DNADB_HOST BUILD
			
			#This could be the admin user, which we don't want for DNADB
			export DNADB_USER=${DNADB_USER:='ensro'}
			#We're not defaulting host, so let's not default port either
			export DNADB_PORT=${DNADB_PORT:='3306'}
			export DNADB_USER=${DNADB_USER:='ensro'}
   			#export DNADB_OPTIONS=' -
			#export DNADB_MYSQLARGS
		fi


		seq_type="${align_type}SEQS"
		found_file=
		seq_path=$(eval "echo \$$seq_type")

		#Test seq path		

		if [ $seq_path ]
		then

			#echo 'seq path defined'
		
			if [ ! -f "$seq_path" ]
			then
				echo "Could not find $seq_type: $seq_path"
			else
				found_file=1
			fi
		fi
	

		#Use default seq file based on validated build?

		if [ ! $found_file ]
		then
			seqs_dir="${ARRAYS_HOME}/${seq_type}"
				
			if [ -d $seqs_dir ]
			then

				if [ ! $found_file ]
				then

					if [ $seq_type = GENOMICSEQS ]
					then
						file_regex=toplevel_*

					elif [ $seq_type = TRANSCRIPTSEQS ]
					then
						file_regex=transcripts.*
					fi


					# Now ask if we want to use any of the files currently present					
					files=()

					for file in $(ls $ARRAYS_HOME/$seq_type/${SPECIES}_${file_regex}.fasta)
					do
					  files=(${files[*]} $file)
					done		
								
					num_files=${#files[@]}
					#This can be an empty string so num_file can be 1	
		  			i=0
				
		
					if [ $files ] && [ $num_files -gt 0 ]
					then
						echo "Would you like to use the following $seq_type file?"
				
                    #Clean REPLY before while just incase we have inherited it from a previous question
						REPLY=

					#Need to evaluate vars here with ""'s as might be undef, which would cause error
	
						while [ "$REPLY" != y  ] && [ $i -lt "$num_files" ]; do
				   		#This causes a hang
						#REPLY=$(AskQuestion "Would you like to use the following $seq_type file? [y|n] $file")
							AskQuestion "${files[$i]} [y|n]"
						
							if [ "$REPLY" = y ]; then
								eval $seq_type=${files[$i]}
								echo ":: Using $seq_type file: ${files[$i]}"
								found_file=1
							fi
						
							i=`expr $i + 1`
						done	
					fi
				fi
			fi

			if [ ! $found_file ]
			then
				AskQuestion "No $seq_type file found, would you like to DumpSeq? [y|n]"

				if [ "$REPLY" = y ]
				then
					DumpSeq $align_type
					found_file=1
				fi
			fi	
		fi

		if [ ! $found_file ]
	   	then
	   		echo "WARNING: Could not set valid $seq_type file, please dump manually or try again"
			exit 1;
     	fi
	done

	echo "GENOMICSEQS:      $GENOMICSEQS
TRANSCRIPTSEQS:   $TRANSCRIPTSEQS

"
#cd $WORK_DIR
	workdir 
#aliases in pipeline.env not picked up?

}


#Could we move this to pipeline.env and make generic?
#Need local mode here

DumpSeq(){
	seq_type=$1
	#should take output_file var here?

	CheckVariables seq_type

	#Validate
	ValidateVariable seq_type VALID_ALIGN_TYPES

	#seqs_var="${seq_type}SEQS"
   	dir="${ARRAYS_HOME}/${seq_type}SEQS"
   	
	#Create SEQS dir with large stripe size

	if [ ! -d $dir ]
	then
		MakeDirs $dir
		lfs setstripe $dir 0 -1 -1
	fi
	
	dump_dnadb_param=

	if [ $DNADB_PASS ]; then
		dnadb_pass_param="--dbpass $DNADB_PASS"
	fi

	jid=


	#These all need to chage to schema build due to 
	#lack of genebuild version and additional patches on the same assembly


	#assm=$(echo $BUILD | sed 's/[a-z]//')


	if [ $seq_type = "GENOMIC" ]; then
		#We don't yet support non-ref seqs here MHC haps etc.
		#Would need to revert to schema_buld as dbname soes not encode patches
		job_name=${SPECIES}_toplevel_${SCHEMA_BUILD}.fasta
		export GENOMICSEQS=${dir}/${job_name}
		echo ":: Dumping softmasked genomic sequence: $GENOMICSEQS"

		#This sometimes runs out of memory on head nodes??
		#Always after chr 6 and before chr 4 for human
		#Even tho mem usage has remained ~ 20%???

		#Do this as CVS co will not be executable
		chmod +x $ANALYSIS_SCRIPTS/sequence_dump.pl
		BackUpFile $GENOMICSEQS
	
 	#Can't do this if we have note made a full copy of the DB!!!
		bsub_cmd="-e ${dir}/genomic_seq.%J.err -o ${dir}/genomic_seq.%J.out -R\"select[mem>3500] rusage[mem=3500]\" -M3500000 "
		cmd="$EFG_PERL $ANALYSIS_SCRIPTS/sequence_dump.pl \
                 -dbhost $DNADB_HOST \
                 -dbuser $DNADB_USER \
                 -dbname $DNADB_NAME \
                 -species $SPECIES \
                 -coord_system_name toplevel \
                 -mask_repeat Dust  \
                 -mask_repeat RepeatMask  \
                 -softmask $MULTI_SPECIES \
                 -filename $GENOMICSEQS $dnadb_pass_param"

	elif  [ $seq_type = "TRANSCRIPT" ]; then
		#CheckVariables now done in _InitEnv



		job_name=${SPECIES}_transcripts.${SCHEMA_BUILD}.fasta
		export TRANSCRIPTSEQS=${dir}/${job_name}
		echo ":: Dumping transcript sequence: $TRANSCRIPTSEQS"

		BackUpFile $TRANSCRIPTSEQS

		#If we are doing a TRANSCRIPT mapping, do we have to define the DNADB if we already have the correct transcript dump?
		#The default DNADB may not be the correct one
		#Need to define DNADB options string if we are doing TRANSCRIPT

		bsub_cmd="-e ${dir}/transcript_seq.%J.err -o ${dir}/transcript_seq.%J.out -R\"select[mem>3500] rusage[mem=3500]\" -M3500000 "




		cmd="$EFG_PERL $EFG_SCRIPTS/export/dump_genes.pl \
				-dbport  $DNADB_PORT \
				-dbuser  $DNADB_USER \
				-dbhost  $DNADB_HOST \
				-dbname  $DNADB_NAME \
				-species $SPECIES \
				$MULTI_SPECIES  \
				-cdna  \
				-stable_id \
				-file $TRANSCRIPTSEQS $dnadb_pass_param"
			
	else
		echo "The seq_type $seq_type is not recognized"
		exit 1;
	fi


	submitJob $job_name "$bsub_cmd" "$cmd"
	
	if [[ "$QUEUE_MANAGER" != Local ]]; then
		echo "Now need to wait or $seq_type dump to finish: jobWait $JOB_ID"
	fi

}


################################################################################
# Func      : RunAlignments
# Desc      : Wrapper method to run the whole genomic and transcript array 
#             mapping and xref pipeline 
# Args [1]  : Optional - Custom arrays.list file
# Args [2]  : Optional - ARRAY_FORMATS for custom arrays.list file
# Return    : none 
# Exception : none
################################################################################

RunAlignments(){
	#Make this take multiple formats?
	#Is this possible to combine with the file?
	#Are we getting confused with logic names here? ProbeAlign ProbeTranscriptAlign?
	
	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	array_formats=
	array_file=
	align_types=
	skip_import=
	usage='usage: RunAlignments [ -f format ]*  [ -t ALIGN_TYPE(e.g. GENOMIC) -l custom_array_list.file  -s(kip import) -h(elp)]'


	#Can we take array_names here too?
	#Is this wise to restrict to arrays within a linked set?


	while getopts ":f:l:t:hs" opt; do
		case $opt in 
	        f  ) array_formats="$array_formats $OPTARG" ;; 
	t  ) align_types="$align_types $OPTARG" ;;
	s  ) skip_import=' -s ' ;;
	l  ) array_file=$OPTARG ;;
	h  ) echo $usage; return 0;;
	\? ) echo $usage; exit 1;;
	esac 
	done


	#Do this here so we don't have to pass to other methods
	SetArrayFormats $array_formats
	
	file_param=
	
	if [[ $array_file ]]; then
		file_param="-l $array_file"
	fi

	#Validate/Set formats and align types
	SetAlignTypes $align_types
	SetArrayFormats $array_formats
	#Do this here so we don't have to pass to other methods
	
	echo "Array Mapping Pipeline Params are:"
    echo "Funcgen DB:      $DB_MYSQL_ARGS"
	echo "Pipeline DB:     $PDB_MYSQL_ARGS"
	echo "DNA DB:          $DNADB_MYSQL_ARGS"
	echo "GENOMICSEQS:     $GENOMICSEQS"
	echo "TRANSCRIPTSEQS:  $TRANSCRIPTSEQS"

	#Don't need to pass array_formats in the following methods
	#Unless we want to change them
	BuildFastas $file_param
	#Build Fastas need updating as it is not fully tested in custom file mode

    SetUpPipeline $skip_import
	
	if [[ ! $skip_import ]]; then
		BackUpTables -t arrays
		ImportArrays 
		ImportWait
	fi
	
    CreateAlignIDs 
    SubmitAlign
	AlignWait
    monitor
	ProbeAlignReport
    
}



################################################################################
# Func      : BuildFastas
# Desc      : Merges array fasta files in to array format specific files for 
#             ImportArrays.  
# Args [1]  : Optional - Custom arrays.list file
# Args [2]  : Optional - array_format for custom arrays.list file
# Return    : none 
# Exception : Exits if custom file specified but not array format specific
################################################################################

BuildFastas(){
    CheckVariables WORK_DIR
	
	echo ":::: BuildFastas $* ::::"
   		
	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	array_list=
	array_formats=

	usage='usage: BuildFormats [ -f format ]*  [ -l custom arrays.list ]'

	while getopts ":f:l:h" opt; do
		case $opt in 
	        f  ) array_formats="$array_formats $OPTARG" ;; 
            l  ) array_list=1 ;;
			h  ) echo $usage; return 0;;
            \? ) echo $usage; exit 1;;
		esac 
	done


	cnt=0
	
	

	#This currently doesn't allow for multiple formats in one array.list
	#We need to either implement format in list file and change this sub to handle
	#Or we need to run spearately for each format, and add format to ARRAY_FORMATS
	#Do we relaly want to have the default used?
	
	#We also need to be able to FlushFasta or DeleteFastas to enable regeneration from scratch?
	#Why don't we just always build them from scratch?
		
	
	
	#echo "array_list is $array_list"
	#echo "array_formats $array_formats"
	
	
	BackUpFile $WORK_DIR/arrays.list
	
	if [ $array_list ]; then
		echo ":: Building $array_format fasta from custom array list: $array_list"

		echo "Not supported in efg yet"
		exit;
		
		if [ ! -f $array_list ]; then
			echo ":: Array list file does not exist: $array_list"
			exit 1;
		fi

		CheckVariables array_formats
   		cat $array_list
		
		format_fasta=$WORK_DIR/arrays.${array_formats}.fasta
		BackUpFile $format_fasta
		ChangeDir $WORK_DIR
		
		#We could do with some healtchecks here
		#AskQuestion if files already present for array_format
		#prompt to clean array and probe tables?
		
		while read path; do

			#Need to split path and format here

			if [[ -f $path ]]; then
				cp $path $WORK_DIR
				file=$(GetFilename $path)
				
				UnzipCatFasta $file $format_fasta
				
				if [ $? -ne 0 ]; then
					cnt=`expr $cnt + 1`
					echo $file >> $WORK_DIR/arrays.list
				fi
			else
				echo ":: Exiting: Could not find $path"
				return 1
			fi
			
		done < $array_list
		
		#Should copy $array_list to WORK_DIR if it isn't the same
	    	    		
    else # ! $array_list
		SetArrayFormats $array_formats
		arraysdir
				
		#We now need to build individual lists for each ARRAY_FORMAT
		#ignoring the TRANSCRIPTSEQS/GENOMICSEQS and maybe README?
		#Restrict ls to only list directories?
		#ls -d? only lists current dir as .?
		#Also account for ARRAY_FORMATS here, which may have been defined in the env.

		if [ ! -z "$ARRAY_FORMATS" ]; then
			echo ":: Building fasta from specified formats: ${ARRAY_FORMATS}"
			list=$ARRAY_FORMATS
		else
			echo ":: Building fasta from default array formats in ${ARRAYS_HOME}"
			list=$(ls)
		fi


		
		for format in $ARRAY_FORMATS; do
			echo ":: Building fasta file for $format"

		#We could do with some healthchecks here
		#AskQuestion if files already present
		#prompt to clean array and probe tables?
		#Or will pipeline fail when we see a duplicate array?
		#We have should have an array_format subdir
			
			format_home=${ARRAYS_HOME}/${format}
			format_fasta=$format_home/arrays.${format}.fasta
			
		#BackUp when we are dealing with the WORK_DIR version
			
		#this is necessary to allow the grep below
			if [[ ! -f ${format_home}/arrays.list ]]; then
				touch ${format_home}/arrays.list
			fi
			
			ChangeDir $format_home
			

			#change this to *fa || *fasta?

			for file in $(ls); do
				
				if [[ $file = *fa ]] || [[ $file = *fasta ]]; then

					if [ ! -d $file ]; then
	
						if [[ $file != arrays.${format}.fasta ]]; then 	#Skip previous arrays.fasta
							
							if [ ! $(grep "^$file" ${format_home}/arrays.list) ]; then
					#This will catch any files which are a sub string of the zipped file
					#So should skip both the zip and the unzipped file
								
					#However this may also mean some arrays maybe missed if they are 
					#sub strings of pre-existing arrays, but only likely if they are added unzipped
					#e.g. first_array.with_suffix.gz
					#first_array.gz
					#This will not cause a problem because of the gz
					#but add it unzipped
					#first_array 
					#Will grep first_array.with_suffix.gz
					#And therefore will be skipped
							
								UnzipCatFasta $file $format_fasta
								
								if [ $? -ne 0 ]; then
									cnt=`expr $cnt + 1`
									#echo $nfile >> $format_home/arrays.list
								#This causes duplicates!
								#We need to build this at the end of build fastas using the individual format lists
								#$format_home/arrays.list is not used anymore?  This is only a record of the last run
								fi
							fi
						fi
					fi
				fi
			done


			#Regenerate softlink in case data/workdir has moved
			if [ -L $WORK_DIR/arrays.${format}.fasta ]; then
				
				rm -f $WORK_DIR/arrays.${format}.fasta
			fi
			
			Execute ln -s $format_fasta $WORK_DIR/arrays.${format}.fasta
			cat $format_home/arrays.list >> $WORK_DIR/arrays.list
			
		    #Change back to species/array home here
			arraysdir
		done
	fi


	#Now just print out some details
	
	if [[ $cnt -eq 0 ]]; then
	    echo ":: No new fasta files found"
	else
		echo ":: Found $cnt new array fasta files"
	fi
	
    cnt=($(wc -l $WORK_DIR/arrays.list))
	cnt=${cnt[0]}
    #echo ":: $cnt array fasta files cat'd to $WORK_DIR/arrays.list:"
	cat $WORK_DIR/arrays.list
	echo ""
	
	workdir
}



################################################################################
# Func      : UnzipCatFasta
# Desc      : Unzips fasta file if required and tests for ^> fasta headers
# Args [1]  : Fasta file path
# Return    : Boolean - true is succesful, false if not fasta file
# Exception : None
################################################################################


#We need to remove this or make it aware of files that are already unzipped
#so we don't get duplication

UnzipCatFasta(){
	
	

    file=$1
	shift 
	format_fasta=$1
	shift

	CheckVariables file format_fasta


	#Need to implement IsCompressed here?
	#Need to grep file for archive or compressed here
	
	#'Zip archive data'
	#'gzip compressed data'
	file_info=$(file $file)
	
    if [[ $file_info = *Zip* ]]; then
		echo "unzipping"
		unzip $file
		rm -f $file
		nfile=$(echo $file | sed 's/\.zip//')

	elif [[ $file = *gzip* ]]; then	
		gunzip $ofile
		nfile=$(echo $file | sed 's/\.gz//')
	else
		nfile=$file
	fi
	


	if [ ! -e $nfile ]; then
		echo "Could not detect unzipped file name for $file, please try BuildFastas again"
		#This should work as we have removed the original file
		exit
	fi
		

    #This awk will not work on large files!!
	#count=$(awk  'BEGIN { RS = "/^>/" }; /^>.*\n^[AGCTUagctu]+\n/' $nfile | wc -l) 
	#count=$(( ($count - 1) /2 ))

	#Cannot do wrapped grep
	#header and non seq lines must match!

	
	dos2unix $nfile 
	count=$(grep -c -E "^>" $nfile)
	#Some Phalanx have N's?!
	non_seq_count=$(grep -c -v -E "^[NAGCTUnacgtu]+$" $nfile)
	


	if [ $count -ne $non_seq_count ]; then
		echo "Headers:$count Non headers:$non_seq_count"
		echo "Found non valid fasta record in $nfile"
		
		grep -E "^>" $nfile > $WORK_DIR/headers.tmp
		grep -v -E "^[NAGCTUnacgtu]+$" $nfile > $WORK_DIR/not_seq.tmp 

		diff $WORK_DIR/headers.tmp $WORK_DIR/not_seq.tmp 

		rm -f $WORK_DIR/headers.tmp $WORK_DIR/not_seq.tmp 

		exit 1
	fi
	
    if [[ $count -gt 0 ]]; then
		echo ":: Found array fasta file: $file"
		echo ":: Contains $count probes"
		fasta=1
		cat $nfile >> $format_fasta
    else
		echo ":: Skipping non-fasta file: $file"
		fasta=0
	fi
	
    return $fasta
}


SetAlignTypes(){
	align_types=$*

	if [[ $align_types ]]; then
		
		for type in $align_types; do
			ValidateVariable type VALID_ALIGN_TYPES
		done

		ALIGN_TYPES=$align_types
		echo ": Setting align types: $ALIGN_TYPES"

	elif [[ ! $ALIGN_TYPES ]]; then
		ALIGN_TYPES=$VALID_ALIGN_TYPES
		echo ": Setting default align types: $ALIGN_TYPES"
	fi
	
	#else we have ALIGN_TYPES which has been previously Set & therefore validated
}



################################################################################
# Func      : SetArrayFormats
# Desc      : Set the ARRAY_FORMATS varible with passed args, or defaults to
#             previously set ARRAY_FORMATS  or the array formats present in 
#             $ARRAYS_HOME   
# Args [*]  : optional - list of valid array formats
# Return    : none 
# Exception : none
################################################################################

SetArrayFormats(){
	formats=$*
	CheckDirs $ARRAYS_HOME
	cwd=$PWD
	arraysdir
	skip=

	
	if [ -z "$formats" ]; then
		#$formats not set

		if [ ! -z "$ARRAY_FORMATS" ]; then
			#$ARRAY_FORMATS previously set
			skip=1
		else
			echo ": Setting default ARRAY_FORMATS in: ${ARRAYS_HOME}"
			tmp=$(ls $ARRAYS_HOME)
			formats=

			for entry in $tmp; do
				
				if [[ ! $entry = GENOMICSEQS ]] && [[ ! $entry = TRANSCRIPTSEQS ]]; then

					if [ -d $entry ]; then
						formats="$formats $entry"
					fi
				fi
			done

			echo ": ARRAY_FORMATS are $formats"
		fi
	fi

	#Validate new formats
	if [ ! $skip ]
	then
	
		for format in $formats; do

		  if [ ! -d $format ]; then
			  echo "You have specified an array format which does not exist in your input dir: $PWD/$format"
			  exit 1
		  fi

		  #we need to capture error here and warn about only having valid format dirs in arrayshoms
		  ValidateVariable format VALID_ARRAY_FORMATS
		done		 

		export ARRAY_FORMATS=$formats
		echo ": Setting ARRAY_FORMATS: $ARRAY_FORMATS"
	fi

	cd $cwd
}


################################################################################
# Func      : SetUpPipeline
# Desc      : Imports pipeline tables, analyses and rules into DB
# Args [1]  : 
# Args [2]  : Optional - ARRAY_FORMAT for custom arrays.list file
# Return    : none 
# Exception : none
################################################################################

#We need to separate the DB generation from this so that we can set up individual analyses whilst we have other stuff running?
#Can we import rules/goals on top of each other?
#We may get a running pipeline picking up the jobs straight away, so we need to make sure evertyhing is in place first. i.e. BuildFastas

SetUpPipeline(){
	echo ":: SetUpPipeline $*"

	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	array_formats=
	skip_import=
	usage='usage: SetUpPipeline [ -f format ]* [ -s(kip import set up) ] [ -h(elp) ]'


	#Need to take align type here too

	while getopts ":f:hs:" opt; do
		case $opt in 
	        f  ) array_formats="$array_formats $OPTARG" ;; 
	s  ) skip_import=' -s ' ;;
	h  ) echo $usage; return 0;;
	\? ) echo $usage; exit 1;;
	esac 
	done



	SetArrayFormats $array_formats

	#Set up analysis and rule config files
	analysis_conf_file=$WORK_DIR/probe_transcript_analysis.conf
	rules_conf_file=$WORK_DIR/probe_transcript_rules.conf
	batch_queue_file=$WORK_DIR/batch_queue.pm

	BackUpFile $analysis_conf_file
	rm -f $analysis_conf_file
	touch $analysis_conf_file

	BackUpFile $rules_conf_file
	rm -f $rules_conf_file
	touch $rules_conf_file


	#Can't add params to analysis here as we don't have access yet, but we should do that somewhere?
	#Could do it in the Runnable DB, but highly redundant, as this will be tested/overwritten for each job i.e. thousands of times.
	

	#Space separate logic names to allow readability in Analysis::Config modules


	#This is not right?
	#We can map AFFY and AFFY_ST together, as they have the same probe design, i.e. 25 mers
	#Just keep separate for now for simplicity
	#We also don't want ProbeTranscriptAlign for Tiling arrays?
	#We need to this dependent on ARRAY_DESIGN?  This is current AFFY or AFFY_ST
	#or CLASS. TILING/EXPRESSION.
	#This concept is a little weird as it's the application of the array, not the array design itlsef, 
	#as you can use TILING designs for WG expression analysis.  Then the idea of probe sets becomes a little
	#awry.

	if [ $warn ]; then
		echo 'TO DO: Need to fix analysis and rules set up dependant on array type, i.e. tiling should not be transcript mapped'
		echo 'NOTE: also need to test for transcript coord_system if we are mapping to transcript?'
	fi

	echo ": Writing analysis and rules config"

    #added null modules to stop warnings on import
	
    #These will not overwrite current entries in analysis, so we may have some old data listed
	#Need to validate this?
	
	#Create empty conf files
	rm -f $rules_conf_file
	touch $rules_conf_file
	rm -f $analysis_conf_file
	touch $analysis_conf_file
	rm -f $batch_queue_file
	touch $batch_queue_file
	
	
	#Set up start of local batch_queue.pm

	echo  "
our @QUEUE_CONFIG = (
" >> $batch_queue_file


	align_condition=

	if [ ! $skip_import ]; then
	   
		echo  "        {
	logic_name => 'Import_Wait',
	batch_size => 1,
	resource   => '',
	retries    => 3,
	sub_args   => '',
	runner     => '',
	cleanup => 'no',         
	},
" >> $batch_queue_file

		align_condition='condition=Import_Wait'

	#Build multi-condition waits first
		echo "[Import_Wait]" >> $rules_conf_file






		for format in $ARRAY_FORMATS; do

			echo  "  	{
	 logic_name => 'Import_${format}_Arrays',
	 batch_size => 1,
	 #These need to be dynamic dependent on number of and type of arrays
	 queue      => \$ENV{'IMPORT_ARRAYS_QUEUE'},
	 resource   => \$ENV{'IMPORT_ARRAYS_RESOURCE'},
	 retries    => 3,
	 cleanup => 'no',         
	},
	
	{
	 logic_name => 'Submit_Import_${format}_Arrays',
	}," >> $batch_queue_file



			echo "condition=Import_${format}_Arrays" >> $rules_conf_file
		done


	#Do we not need to add another rule here for the Submit jobs to wait for Import_Wait?
	#Or is this what the Import_Wait accumulator is doing?
	
	
		echo "[Import_Wait]
module=Accumulator
input_id_type=ACCUMULATOR
">> $analysis_conf_file
	
	else
		echo ": Skipping Import setup"
	fi

	
	echo "    {
	 logic_name => 'ProbeAlign_Wait',
	 batch_size => 1,
	 resource   => '',
	 retries    => 3,
	 sub_args   => '',
	 runner     => '',
	 cleanup => 'no',         
	},
    {
	 logic_name => 'Align_Wait',
	 batch_size => 1,
	 resource   => '',
	 retries    => 3,
	 sub_args   => '',
	 runner     => '',
	 cleanup => 'no',         
	},"   	>> $batch_queue_file


	#Set ucfirst align types
	align_types=$(echo $ALIGN_TYPES | tr [A-Z] [a-z] | perl -ane ' foreach $wrd ( @F ) { print ucfirst($wrd)." "; } print "\n" ; ')
	#Default to ImportWait
	talign_condition=$align_condition

	

	#We could actually split this our such that we can
	#deal with each format independantly
	#Hence we can finish formats individually and progress to 
	#ProbeAlignReport and RunTranscriptXrefs instead of having to wait for bigger formats

	if [[ $align_types = *Genomic* ]] && [[ $align_types = *Transcript* ]]; then
		talign_condition='condition=ProbeAlign_Wait'
		echo "[ProbeAlign_Wait]
condition=${format}_ProbeAlign" >> $rules_conf_file	

		echo "[ProbeAlign_Wait]
module=Accumulator
input_id_type=ACCUMULATOR
">> $analysis_conf_file

		#probe_align_wait=1
	fi
	


	#Set resource requirements to enable DB throttling
	
	resource="-R \"select[(${DB_HOST_LSFNAME}<=700)"
	
	#is duration here 
	#each individual job takes around 30 mins
	#to run locally
	#but these are batched jobs, 5 in series?
	#however these seem to run much faster on the farm nodes
	#taking 20 mins for all 5!

	rusage="rusage[${DB_HOST_LSFNAME}=12:duration=10"
	#measured load for this job will only ever be 13
	#10 for active query and 1 each connection to pipeline, core or efg DB.

	if [[ $DB_HOST_LSFNAME != $DNADB_HOST_LSFNAME ]]; then
		resource="${resource}&&(${DNADB_HOST_LSFNAME}<=700)"
		rusage="$rusage,${DB_HOST_LSFNAME}=13:duration=10]"
	fi

	resource="${resource}] ${rusage}]\""



	for format in $ARRAY_FORMATS; do
	
		if [[ $align_types = *Genomic* ]]; then
			echo  "  
	{
	 logic_name => 'Submit_${format}_ProbeAlign',
	},
		
	{
	 logic_name => '${format}_ProbeAlign',
	 batch_size => 5,
	 retries    => 3,
	 #Set temporarily as a few jobs were failing on memory limits
     #also set batch_size to 1 if you are only rerunning a handful of job
	 #Implement ${format}_ProbeAlign_RESOURCE/ALIGN_MEM_FACTOR?
	 #resource   => '-R\"select[mem>6000] rusage[mem=6000]\" -M 6000000',
     resource    => '$resource',
	 cleanup => 'no',
	},"  >> $batch_queue_file
			
		fi


		if [[ $align_types = *Transcript* ]]; then
			echo "
	{
	 logic_name => 'Submit_${format}_ProbeTranscriptAlign',
	},
				   
	{
	 logic_name => '${format}_ProbeTranscriptAlign',
	 batch_size => 7,#No point in reducing this as we'll be waiting for the ProbeAlign jobs
	 retries    => 3,
	 resource    => '$resource',
	 cleanup => 'no',         
	}," >> $batch_queue_file

			#echo "condition=${format}_ProbeTranscriptAlign" >> $rules_conf_file	

		fi
	done


	#Now finish of batch_queue.pm	
	echo -e ");\n1;" >> $batch_queue_file


	#Now add general AlignWait
	#dependant on last alignment analysis
	#Can now either be Genomic or Transcript
	
	echo "[Align_Wait]" >> $rules_conf_file
	atype=''
	
	if [[ $align_types = *Transcript* ]]; then
		atype='Transcript'
	fi


	for format in $ARRAY_FORMATS; do
		echo "condition=${format}_Probe${atype}Align" >> $rules_conf_file	
	done

		

	echo "[Align_Wait]
module=Accumulator
input_id_type=ACCUMULATOR
">> $analysis_conf_file


	#Need to account for align type here!
	#We are currently writing analyses for both
	
	for format in $ARRAY_FORMATS; do

		if [ ! $skip_import ]; then
		
			echo "[Import_${format}_Arrays]
module=ImportArrays
input_id_type=PROBE_SET

[Submit_Import_${format}_Arrays]
input_id_type=PROBE_SET
" >> $analysis_conf_file
			
			echo "
[Import_${format}_Arrays]
condition=Submit_Import_${format}_Arrays
" >> $rules_conf_file
		fi
		

		if [[ $align_types = *Genomic* ]]; then

		echo "[Submit_${format}_ProbeAlign]
input_id_type=PROBE_CHUNK

[${format}_ProbeAlign]
program=exonerate
program_version=$EXONERATE_VERSION
program_file=$EXONERATE_PATH
module=ProbeAlign
input_id_type=PROBE_CHUNK
" >> $analysis_conf_file

		
		echo "
[${format}_ProbeAlign]
$align_condition
condition=Submit_${format}_ProbeAlign
" >> $rules_conf_file

		fi

		if [[ $align_types = *Transcript* ]]; then

			echo "[Submit_${format}_ProbeTranscriptAlign]
input_id_type=PROBE_TRANSCRIPT_CHUNK

[${format}_ProbeTranscriptAlign]
program=exonerate
program_version=$EXONERATE_VERSION
program_file=$EXONERATE_PATH
module=ProbeAlign
input_id_type=PROBE_TRANSCRIPT_CHUNK
" >> $analysis_conf_file

			echo "[${format}_ProbeTranscriptAlign]
$talign_condition
condition=Submit_${format}_ProbeTranscriptAlign
" >> $rules_conf_file

		fi 

		
	done



	CreatePipelineTables
   


	#Could do with testing for files here
	echo ": Importing analysis config" 

	#EFG
    Execute $EFG_PERL $PIPELINE_SCRIPTS/analysis_setup.pl $PDB_SCRIPT_ARGS -read -file $WORK_DIR/probe_transcript_analysis.conf 

    echo ": Importing rules config" 
	#EFG
    Execute $EFG_PERL $PIPELINE_SCRIPTS/rule_setup.pl  $PDB_SCRIPT_ARGS -read -file $WORK_DIR/probe_transcript_rules.conf

	#Could we test output for 'Not storing' This will not detect whether there are other rules in the DB

	if [ $warn ]; then
		echo "NOTE: Need to clean analysis table? Or is this a setup problem"
		echo "NOTE: check for genomic seq here, are all the chrs present in $GENOMICSEQS"
		#This would requite getting all toplevel (inc non-ref) and greping file
		#Is this another case for PipeLineHelper.pm
	fi
	
    Execute $EFG_PERL $PIPELINE_SCRIPTS/setup_batchqueue_outputdir.pl 
	
	#This will always fail and return 0?!
	#As we don't have any input_id written for the accumulator yet
	#This needs to be done after we have created the input ids for ImportArrays

	echo "NOTE: The following 'analysis Submit_*' and accumulator warnings can largely be ignored" 

	CheckPipelineSanity

	echo ""

	}






################################################################################
# Func      : ImportArrays
# Desc      : Collapses arrays of related formats in to unique probe records 
#             based on probeset and sequence identity
# Params    : -f optional - Array format e.g. AFFY_UTR, can be specified more than once
#             -r rollback - Rollback flag remove IMPORTED status from each array to enable
#             rollback in ImportArrays.
# Return    : none 
# Exception : Exits if ARRAY_FORMATS not defined
################################################################################

#Need to add no collapse/straight import functionality for nr_arrays i.e. ST arrays

ImportArrays(){
	echo ":::: ImportArrays $* ::::"

	formats=
	rollback=

	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1
	usage='usage: ImportArrays [ -f format ]*  [ -r(ollback) ] [ -h(elp) ]'

	while getopts ":f:rh" opt; do
		case $opt in 
	        f  ) formats="$formats $OPTARG" ;; 
            r  ) rollback=1 ;;
			h  ) echo $usage; return 0;;
            \? ) echo $usage; exit 1;;
		esac 
	done

			
	SetArrayFormats $formats
	CheckVariables ARRAY_FORMATS ARRAYS_HOME
			
	echo "Need to check pipeline sanity here?"
	echo "Need to set up separate pipeline DB"
	echo "Should clean up job_status here too?"

	#pipeline_sanity.pl -dbhost $DB_HOST -dbport $DB_PORT -dbuser $DB_USER -dbpass $DB_PASS -dbname $DB_NAME
	#This is currently barfing due to presence of other analyses which do not have modules.
	
	workdir
	
	#Now done in wrapper
	#BackUpTables arrays
	
	#We need to write the config file for each format here.
	cnt=0

	#Create input_ids and config file for ImportArrays
	#Let's have 1:genome etc to maintain some continuity
	#export ARRAY_FORMAT_FILE=$WORK_DIR/ImportArrays.config
	#BackUpFile $ARRAY_FORMAT_FILE
	#rm -f $ARRAY_FORMAT_FILE
	#touch $ARRAY_FORMAT_FILE


	#Need to sub relevant parts of this loop to enable TestImportArray?
	#Would not write or bsub, but then were not testing storage code
	#This loop only works as the rule manager finishes after every submit?
	#Or does it...is this running in series?
	#If so we need to kick off with no -analysis flag 
	#And get it to stop after Import_Wait somehow, so we can create the Align IDs etc
	#Or should we just put this in a RunnableDB and then pipeline the whole thing?

	for format in $ARRAY_FORMATS; do
		logic_name="Import_${format}_Arrays"

		#How does this not result in duplicate imports?
		CleanJobs -l $logic_name
		
		if [[ $rollback ]]; then
			
			SetArrayNames $format
			
			if [[ $ARRAY_NAMES ]]; then
				
				names=$(echo $ARRAY_NAMES | sed "s/ /', '/")	
				
				sql="DELETE s from status s, status_name sn, array_chip ac where ac.name in('${names}') and ac.array_chip_id=s.table_id and s.table_name='array_chip' and s.status_name_id=sn.status_name_id and sn.name='IMPORTED'";
				echo ": Deleting IMPORTED status for: $ARRAY_NAMES"
				Execute echo $sql | mysqlefg
			fi
		fi
		
		
	    #Do we want to do this for every run, what if we've had a failure and we're just re-running
	    cmdline="make_input_ids $PDB_SCRIPT_ARGS -single -logic_name Submit_Import_${format}_Arrays -single_name ${format}:genome"
		echo ": Creating Import_${format}_Arrays input IDs" 
		Execute $cmdline
		export RAW_FASTA=$WORK_DIR/arrays.${format}.fasta
	
		#Set musage
		SetMusage -a "$ARRAYS_HOME/$format/arrays.${format}.fasta" -f $IMPORT_MEM_FACTOR -b $IMPORT_MEM_BASE

		musage_k=$(($musage * 1000))
	   
		if [ $musage -gt $MAX_HUGE_MEM ]; then
			echo "Expected memory usage exceeds that available on $HUGEMEM_QUEUE queue: $musage. Modify rsuage generation?"
			exit 1;
		elif [ $musage -gt $MAX_NORMAL_MEM ]; then
			queue=$HUGEMEM_QUEUE
			echo "Need to move things to home dir for $HUGEMEM_QUEUE jobs. This will go away when we import ST arrays directly?"
			exit 1;
		else
			queue=$NORMAL_QUEUE
		fi

		
		#This has to be double quotes or it get's ignored!!???
		#Weird behaviour with -R options, can sometimes also truncate see RunTranscriptXrefs
		resource="-R \"select[mem>${musage}] rusage[mem=${musage}]\" -M $musage_k"
		
		#We need to copy input files from lustre to somewhere where hugemem can see it
		#Also change respective vars and handle moving data back to lustre after ImportArrays
		#Or shall we move everything for Import and then move everything back later?
		
		
		export IMPORT_ARRAYS_QUEUE=$queue
		export IMPORT_ARRAYS_RESOURCE=$resource
		echo ":: Import_${format}_Arrays bsub params: -q $IMPORT_ARRAYS_QUEUE $IMPORT_ARRAYS_RESOURCE"
		
        #We are just running the pipeline through one cycle to submit here
		#This is because we need to create the input_ids for each ProbeAlign step before we carry on
		#We should create another RunnableDB to do this automatically so we can run the pipeline all in one go!
		cmdline="$EFG_PERL $PIPELINE_SCRIPTS/rulemanager.pl $PDB_SCRIPT_ARGS -once -analysis $logic_name -input_id_type PROBE_SET"
		#echo ":: Import_${format}_Arrays: $cmdline"


		#Archive the names files before Executing
		BackUpFile "${WORK_DIR}/arrays.${format}.names"
		Execute $cmdline
		
		
		#THis is running with jobname like nj_arrays_mus_musculus_funcgen_51_37d:default. How do we replace the default with the analysis? Look in rulemanager. This should have used the analysis!!		
	done
		
	
	echo "Completed incremental submission of ImportArrays for formats: $ARRAY_FORMATS"
		
}



################################################################################
# Func      : ImportWait
# Desc      : Accumulator step to wait for all ImportArrays jobs
# Arg[1..n] : None
# Return    : None 
# Exception : None
################################################################################


ImportWait(){
	echo ":: ImportWait $*"
	#add opts here for -force_accumulators and -once?

	cmdline="$EFG_PERL $PIPELINE_SCRIPTS/rulemanager.pl $PDB_SCRIPT_ARGS -force_accumulators -analysis Import_Wait"

	#echo ":: ImportWait $cmdline"

	Execute $cmdline
	
	#echo "Can we cat the AFFY and AFFY_ST arrays_nr files?"
	#As these should use the same rules?  This would mean running with one logic name for both
	#Just keep separate for now

	#Now cat the array names files output by ImportArrays
	#This is to enable any ProbeFeature/IDXref/UnmappedObject deletes

	#export NAMES_FILE="${WORK_DIR}/arrays.names"

	#ArchiveFile NAMES_FILE

	#CheckVariables ARRAY_FORMATS

	#for format in $ARRAY_FORMATS; do
	#	format_names="${WORK_DIR}/arrays.${format}.names"
	#	cat $format_names > $NAMES_FILE
	#done
		



	#We could validate array names or at least number of arrays here?
	echo "Finished waiting...ready for CreateAlignIds [ formats ]"
}



SetArrayNames(){
	formats=$*

	format_files=

	if [ ! -z "$formats" ]; then

		for format in $formats; do
			format_files="${WORK_DIR}/arrays.${format}.names $format_files"
		done
	else
		format_files="${WORK_DIR}/arrays.names"
	fi


	#need to change this so it uses the format files instead of this cat'd file?
		#when do we ever use the cat'd file?

	ARRAY_NAMES=

	for file in $format_files; do

		if [ ! -f $file ]; then
			echo -e "Could not find ARRAY_NAMES file:\t$file"
			echo -e 'You either need to:\n\tGenerateArrayNames from the DB\n\tMigrate from a previous workdir\n\tor ImportArrays'
		else
			
			while read name; do
				ARRAY_NAMES="$name $ARRAY_NAMES"			
			done < $file
		fi
	done


	#Can't do this as we may want to call this format by format
	#and some may not have names files yet?
	#i.e. When ImportArrays has failed
	#if [ -z $ARRAY_NAMES ]; then
	#	echo "Found no valid array names"
	#	exit 1;
	#fi

	#So must test return value in caller		

}



################################################################################
# Func      : CreateAlignIDs
# Desc      : Creates ProbeAlign or ProbeTranscriptAlign input IDs
# Opt -t    : Mandatory - Align types e.g. GENOMIC and/or TRANSCRIPT
# Opt -f    : Optional - Array formats e.g. 'AFFY AFFY_ST', defaults will be used if absent
# Opt -c    : Optional - Chunks, chunks will be calculated if omitted
# Arg[3]    : would need to parse opts
#			  ARRAY_DESIGN????? Would need to create different sets of IDs for different designs
#             This would probably screw up  the chunking based on the IDs, would also need different source file
#             Hence another config file would need to be written. 
# Return    : none 
# Exception : None
################################################################################

CreateAlignIDs(){
	
	echo ":::: CreateAlignIDs $*"

	align_types=
	formats=
	chunk_arg=
	usage='usage: CreateAlignIDs [ -t align_type ]* [ -f format ]* [-c chunks] [ -h(elp) ]'

	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	while getopts ":t:f:c:h" opt; do
		case $opt in 
			t  ) align_types="$align_types $OPTARG" ;;
		f  ) formats="$formats $OPTARG" ;; 
		c  ) chunk_arg=$OPTARG ;;
		h  ) echo $usage; return 0;;
		\? ) echo $usage; exit 1;;
		esac 
	done

	#No shift required as we don't expect any other args
    #shift $(($OPTIND - 1))

	SetAlignTypes $align_types
	SetArrayFormats $formats
	#We need to back up the input ids before deleting and recreating
	



	for format in $ARRAY_FORMATS; do
	
		#If we are re-initialising the env
		#and have only built the fastas using a subset of
		#of formats, this will currently try and create IDs for all available formats
		#and will give errors as the fasta will not be found
				#Just do the fasta count once if required
		

		NR_FASTA="$WORK_DIR/arrays_nr.${format}.fasta"
		
		if [[ ! -f $NR_FASTA ]]; then
			echo "Skipping $format, NR_FASTA not found: $NR_FASTA"
		else
			
			if [ ! $chunks_arg ]; then
				num_probes=`expr $(wc -l $NR_FASTA | sed "s| $NR_FASTA||") / 2`
				echo ": Found $num_probes $format NR probe sequences"
			fi
			
			for align_type in $ALIGN_TYPES; do	 
				
  				if [ $align_type = GENOMIC ]; then
					id_type='PROBE_CHUNK' 
					logic_name="${format}_ProbeAlign"
					
				elif [ $align_type = TRANSCRIPT ]; then
					id_type='PROBE_TRANSCRIPT_CHUNK'
					logic_name="${format}_ProbeTranscriptAlign"
				fi
				
            #Delete the old input IDs here before importing!?
				CleanJobs -l $logic_name

            #Calculate the number of chunks
				if [ $chunk_arg ]; then
					chunks=$chunk_arg
				else
                #Tune the chunk numbers wrt array type and align type
					
                #This is roughly based on nimblegen probes being 50-60mers
	    	    #And affy probes being roughly 25mers				
					if [[ $ARRAY_DESIGN = NIMBLEGEN_TILING ]]; then
						probes_per_chunk=8000
					    #May need to change this
					else
					#probes_per_chunk=8000
					#Takes ~ 7 mins for Genomic ProbeAlign
					#On head node
					#probes_per_chunk=30000
					#With batch size 5, this should take between 2.5-5 hours!!
					#This is taking >3 hours?
					#On a fully loaded node
					#probes_per_chunk=10000
					#Takes >1.5 hours on fully loaded node
					#probes_per_chunk=5000

						probes_per_chunk=3000
					fi

					if [[ $logic_name = *ProbeTranscriptAlign ]]; then	
					#Multiply for transcripts, 10 for now
						probes_per_chunk=`expr $probes_per_chunk \* 10`
					fi
	
					#Add 1 to avoid running with 0
					chunks=$((($num_probes / $probes_per_chunk) + 1))	
					
				fi

		
				echo ": Creating $chunks $logic_name IDs with $probes_per_chunk probes/chunk (e.g. 1:$chunks:$format)"	
				anal_id=$(QueryVal PIPELINE "select analysis_id from analysis where logic_name='Submit_${logic_name}'")
				cmd="for \$i(1..$chunks){ print \"insert into input_id_analysis(input_id,analysis_id,input_id_type) values (\\\"\$i:$chunks:$format\\\",$anal_id,\\\"$id_type\\\"); \n\" }"	
			
			    #Execute "perl -e '$cmd' > ${WORK_DIR}/probe_input_ids.sql"
			    #Cannot get Execute to work for this perl -e
				perl -e "$cmd" > ${WORK_DIR}/probe_input_ids.sql
			    #Execute "mysql $PDB_MYSQL_ARGS < ${WORK_DIR}/probe_input_ids.sql"
				mysqlpipe < ${WORK_DIR}/probe_input_ids.sql

			#Catch error here
			#This is currently causing problems as you can't have the same input_id for two different analyses!!!
            #We now need to parse the input_id differently in ProbeAlign, or is this done in ExonerateProbe?
	done
	fi
	done
	
	echo ": Ready for SubmitAlign"
	
	}


#Do all resources in one go 
#So this can be done by CreateAlignIDS
#And then we can have separate align waits
#and can therefore kick of probetranscriptalign RunnableDB asynchronously
#if we restart the env after failure?
#Therefore should be done in SubmitAlign(which well eventually just run straight through to probe2transcript)


SetResources(){
	echo ": SetResources $*"
	#We may not have Transcript seqs as we may only be doing genomic mapping?
	#Do for all now and change this if we want to include genomic only mapping
	#We have other thing we need to fix for this?

	types=
	formats=
	chunk_arg=
	usage='usage: SetResources [ -t resource_type  e.g. TRANSCRIPT, GENOMIC or XREF ]* [ -f format ]* [ -h(elp) ]'

	

	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	while getopts ":t:f:c:h" opt; do
		case $opt in 
			t  ) types="$types $OPTARG" ;;
		f  ) formats="$formats $OPTARG" ;; 
		h  ) echo $usage; return 0;;
		\? ) echo $usage; exit 1;;
		esac 
	done

	SetArrayFormats $formats
	types=${type:="XREF $VALID_ALIGN_TYPES"}

	#Now we need multipliers $XREF_RESOURCE_MULTIPLIER $TRANSCRIPT_RESOURCE_MULTIPLIER $GENOMIC_RESOURCE_MULTIPLIER
	

	for type in $types; do
		

		#We want to set eg 
		#GENOMIC_AFFY_UTR_MEM_MB
		#GENOMIC_AFFY_UTR_MEM_K
		seq_type=

		if [ $type = XREF ]; then
			seq_type=TRANSCRIPT
		fi

		seq_type=${seq_type:=$type}
		seq_file=$(eval "echo \$${seq_type}SEQS")

		#We need to get number of chunks here


		#We need to set a var name using a var here using eval
		
		eval ${type}_${format}_MEM_K=$mem

		eval ${type}_${format}_MEM_MB=$(($mem / 1000 | bc))


	done
}


#Add verbose option here?

TestProbeAlign(){
	echo ":: TestProbeAlign $*"

	align_type=
	formats=
	chunk=
	write=
	logic_name=
	usage='usage: TestProbeAlign -t GENOMIC|TRANSCRIPT -f "array_format" -c "chunk" [ -w(rite flag) ] [ -h(elp) ]'

	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	while getopts ":t:f:c:wh" opt; do
		case $opt in 
			t  ) align_type=$OPTARG ;;
	        f  ) format=$OPTARG ;;
            c  ) chunk=$OPTARG ;;
			w  ) write=' -write ' ;;
			h  ) echo $usage; return 0 ;;
		    \? ) echo $usage; return 1 ;;
		esac 
	done



	CheckVariablesOrUsage "$usage" format align_type	
	ValidateVariableOrUsage "$usage" format VALID_ARRAY_FORMATS
	ValidateVariableOrUsage "$usage" align_type VALID_ALIGN_TYPES
	
	
	#Check nr fasta file
	CheckFile "${WORK_DIR}/arrays_nr.${format}.fasta"
	
	
	#Set the logic name
	if [[ $align_type = "GENOMIC" ]]; then
		CheckVariables GENOMICSEQS
		logic_name="${format}_ProbeAlign";
	else
		CheckVariables TRANSCRIPTSEQS
		logic_name="${format}_ProbeTranscriptAlign";	
	fi
	
	#Get the chunk count for this analysis
	chunks=$(QueryVal PIPELINE "select count(*) from input_id_analysis i, analysis a where a.logic_name='Submit_${logic_name}' and a.analysis_id=i.analysis_id")
	
	if [ $chunks = 0 ]; then
		echo "Could not find any input_ids for analysis $logic_name";
		exit 1;
	fi
	
	#Take the first chunk if it is not defined
	chunk=${chunk:=1}
	input_id="${chunk}:${chunks}:$format"
	
	#Need to query for inpud_id 
	
	cmd="$ANALYSIS_SCRIPTS/test_RunnableDB $PDB_SCRIPT_ARGS -logic_name $logic_name -input_id $input_id $write"
	echo $cmd
	
	time $EFG_PERL $cmd	
}

#This should only be used to test the parsing of the fasta file
#-write should not be specific 
#Altho this will still write Array and ArrayChip info to the DB
#This will automatically be rolled back if not fully 'IMPORTED'


TestImportArrays(){
	echo "TestImportArrays $*"

	format=
	write=

	#Need to test for jsut one format
	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	usage='usage: TestImportArrays -f format [ -w(rite) ] [ -h(elp) ]';

	while getopts ":f:wh" opt; do
		case $opt in 
	        f  ) format=$OPTARG ;; 
w  ) write=' -write ';;
h  ) echo $usage; return 0;;
\? ) echo $usage; return 1;;
		esac 
	done

	
	CheckVariablesOrUsage "$usage" format
	#Will this still echo the warning from ValidateVariable
	ValidateVariableOrUsage "$usage" format VALID_ARRAY_FORMATS

	
	#No write here!
	#We should turn on debugging here to test the output
	#-verbose only works in the test script not in the RunnableDB!
	#Probably best not to over load the Runnable with iterative tests for verbose

	cmd="$ANALYSIS_SCRIPTS/test_RunnableDB $PDB_SCRIPT_ARGS -logic_name Import_${format}_Arrays -input_id ${format}:genome $write"
	echo $cmd

	time $EFG_PERL $cmd

}


CleanImportRules(){
	echo ": CleanImportRules $*"

	sql='delete rc from rule_conditions rc where rc.rule_condition like "Import_%"'
	echo $sql | mysqlpipe
	sql='delete rg from rule_goal rg, analysis a where rg.goal=a.analysis_id and a.logic_name like "Import_%"'
	echo $sql | mysqlpipe
	#Do we need to add input_id_analysis and input_is_type_analysis here?
}

#Do we need ReSubmitAlign?
#This would set batch size to user defined amount
#



#Should this also CleanJobs?
#Can we not do this from the RunnableDB and use the roll_back_ArrayChip method?


RollbackArrays(){
	echo ":: RollbackArrays $*"

	formats=
	keep_xrefs=
	mode=
	vendor=

	#Need to test for jsut one format
	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	#We could use the Pipeline stages here Instead? e.g. Import, Align, TranscriptXrefs? 

	usage='usage: RollbackArrays [-f format]+ [ -m mode e.g. probe2transcript | ProbeAlign | ProbeTranscriptAlign | probe_feature | probe(default)  -d(force delete) -k(eep probe2transcript xrefs) -v (vendor, only required if format exists on more than vendor platform, specify one format only) -h(elp) ]';

	while getopts ":f:v:m:adkh" opt; do
		case $opt in 
	        f  ) formats="$formats $OPTARG" ;; 
            v  ) vendor=$OPTORG ;;
            m  ) mode=" -m $OPTARG " ;;
            d  ) force=' -force ' ;;
            k  ) keep_xrefs=' -keep_xrefs ' ;;
            h  ) echo $usage; return 0;;
            \? ) echo $usage; return 1;;
		esac 
	done

	array_names=


	if [ "$vendor" ]; then
		vendor="-v $vendor"
	fi

	for format in $formats; do
		
		ValidateVariable format VALID_ARRAY_FORMATS
		SetArrayNames $format
		

	#Try and get the vendor from the DB
	#if [ ! $vendor ]; then
	#	vendor=$(_getVendorByFormat $format)
	#	if [ $? -ne 0 ]; then
	#		echo -e $vendor
	#		return 1
	#	fi
	#fi
	#vendor no longer mandatory

		
		if [[ ! $ARRAY_NAMES ]]; then
			echo -e "Please generate:\t${WORK_DIR}/arrays.${format}.names"

		#What do we want to do here?
		#Allow delete all in DB for given format
		#Or do we want to force migration of arrays.FORMAT.names files
		#To ensure that we only Rollback those which we have dealt with in this environment?

		#Can we add a MigrateImportFiles
		#This would include all names files + ?
		#There is no guarantee these would be full and complete
		#So best to manually curate?
		#Or can we have GenerateArrayNames?
		#Which will grab these from the DB?
		#These will be overwritten if we reimport some new arrays

			return 1
		fi

		array_names="$array_names $ARRAY_NAMES"
		
	#We are not passing any array names here!
	#Which means -v is being gobbled up
	#Did we forget to check something in here?
	#No! If we use a new DB here, then there will be no array.names available
	#We need a migrate function here!!

	#Or maybe we just need to set this optionally and get arrays based on format only?
	done

	Execute $EFG_PERL $EFG_SRC/scripts/rollback/rollback_array.pl $DB_SCRIPT_ARGS $DNADB_SCRIPT_ARGS -s $SPECIES -a $array_names $vendor $mode $force $keep_xrefs
}


GenerateArrayNames(){

	echo ":: GenerateArrayNames $*"

	formats=
	delete_file=
	all=

	#This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	#We could use the Pipeline stages here Instead? e.g. Import, Align, TranscriptXrefs? 

	usage='Description:\tGenerates a file containing a list of array names for a given format.  This list is required for RollbackArrays.\nUsage:\t\tGenerateArrayNames -f format [ -a(ll formats) -d(elete existing file) -h(elp) ]';

	while getopts ":f:adh" opt; do
		case $opt in 
	        f  ) formats=$OPTARG ;; 
            a  ) all=1 ;;
            d  ) delete_file=1 ;;
            h  ) echo -e $usage; return 0;;
			\? ) echo -e $usage; return 1;;
		esac 
	done


	if [ ! $all ] && [ ! $formats ]; then
		echo "You must supply a -f(ormat) or -a(ll) formats option"
		return 1

	elif [ $all ] && [ $formats ]; then
		echo "You must supply either a -f(ormat) or -a(all) formats option"
		return 1
	elif [ $all ]; then
		SetArrayFormats
		formats=$ARRAY_FORMATS
	fi


	
	for format in $formats; do
		
		ValidateVariable format VALID_ARRAY_FORMATS
		
		format_file=${WORK_DIR}/arrays.${format}.names 
		
		ARRAY_NAMES=$(QueryVal OUT "select name from array where class='$format'")
		echo ":"
		echo ": ARRAY_NAMES for $format are $ARRAY_NAMES"
		

		if [ -z "$ARRAY_NAMES" ]; then
			echo ": WARNING: Could not find any ARRAY_NAMES for $format"
		fi

	#Need to test if we are overwriting file here
		
		if [[ -f $format_file ]]; then
			
			if [ ! $delete_file ]; then
				echo -e "Found existing file:\t$format_file\nSpecify -d to overwrite"
				return 1
			else
				rm -f $format_file
			fi
		fi
		
		echo ": ARRAY_NAMES file: $format_file"

		for name in $ARRAY_NAMES; do
			echo $name >> $format_file
		done
	done

	echo ": Please edit files if required"

}


SubmitAlign(){	
	echo ":: SubmitAlign $*"

	OPTIND=1

	#We could use the Pipeline stages here Instead? e.g. Import, Align, TranscriptXrefs? 

	usage='usage: SubmitAlign [ -h(elp) ]\
Simply submits the ProbeAlign jobs defined by CreateAlignIDs';

	#We could add logic_name/analysis here.

	while getopts ":m:h" opt; do
		case $opt in 
			h  ) echo $usage; return 0;;
            \? ) echo $usage; return 1;;
		esac 
	done


	echo "Checking for DumpSeq jobs and external_dbs"
	assm=$(echo $BUILD | sed 's/[a-z]//')

	for atype in $ALIGN_TYPES; do
		

		if [[ $atype = TRANSCRIPT ]]; then
			atype=Transcript
			job_name="${SPECIES}_transcripts.${BUILD}.fasta"
		else
			atype=Genome
			job_name="${SPECIES}_toplevel_${assm}.fasta"
		fi

				
		edb_id=$(QueryVal OUT "SELECT external_db_id from external_db where db_name='${SPECIES}_core_${atype}' and db_release='${SCHEMA_BUILD}'")
	
		if [ ! $edb_id ]; then
			echo -e "Inserting external_db:\t${SPECIES}_core_${atype}\t${SCHEMA_BUILD}"
			echo "INSERT into external_db (db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) values('${SPECIES}_core_${atype}', '${SCHEMA_BUILD}', 'KNOWNXREF', 1, 5, 'Ensembl${atype}', 'MISC')" | mysqlefg
		fi

		checkJob $job_name
		#Sets JOB_ID
		
		if [ $JOB_ID ]; then
			echo -e "You still have a DumpSeq job running:\t$job_name"
			echo -e "Waiting for job($JOB_ID)"
			jobWait $JOB_ID
			
			if [ $? -ne 0 ]; then
				echo -e "$job_name failed with exit code $?"
				exit;
			fi
		fi
	done




	#This will kick off the rest of the pipeline
	#So if we have cleaned the import job details
	#Then it will try and do the collapse step again.
	#Use CleanImportRules here just to be safe?
	CleanImportRules

	#sql to remove import rules
	#delete rc from rule_conditions rc where rc.rule_condition like "Import_%";
	#delete rg from rule_goal rg, analysis a where rg.goal=a.analysis_id and a.logic_name like "Import_%";

	#CleanArrays here? No need, will fail if they have already been imported 
	#or will be rolled back if not(and fail if associated xrefs found)

		#Need to add DB throttling here too
		#rusage[mem=${musage}:$DBHOST=10:duration=10:decay=0]
	#rusage="mem=${musage}"

		#if [ $RUSAGE_DB_HOST ]; then
		#	rusage="${rusage}:"
		#fi



	Execute $EFG_PERL $PIPELINE_SCRIPTS/rulemanager.pl $PDB_SCRIPT_ARGS

	#This is giving weird batching behaviour where the stdout files
	#for some jobs list the analysis of another
	#This is due to the BatchQueue config not having logic name config for each format.


	return;

	#We can only run this once, so we need to remove this loop and check conflicting options
	


	align_types=
	formats=
	chunk_arg=

	echo "before getopts $*"


	#If we reource the env whilst already in the env, then this fails to set vars and align_types check fails!?

	#getopts does not run if we resource the environemnt whilst we are still in the environment!!!

	while getopts ":t:f:c:" opt; do
		case $opt in 
			t  ) echo "found $opt with $OPTARG"; align_types="$align_types $OPTARG" ;; #align_types=$(GetOptArgArray $*) ;;
	        f  ) formats=$OPTARG ;;
            c  ) chunk_arg=$OPTARG ;;
			? ) echo 'usage: CreateAlignIDs -t "Align types" [-f "array formats"] [-c chunks]'; exit 1;;
		esac 
	done


	echo "vars are types: $align_types formats: $formats chunks: $chunk_arg"

	CheckVariables align_types
	SetArrayFormats $formats

	echo ":: Align types: $align_types"



    #Could do with knowing probe size for these calculations
    #Can we pass this as an arg?

	#We should delete the input IDs here before importing!
	

	echo "Submitting ProbeAlign jobs"

	for format in $ARRAY_FORMATS; do
	
		#Need to set this in the runnable as we have to submit in one go.
		#Same with TARGET_SEQS?
		NR_FASTA="$WORK_DIR/arrays_nr.${format}.fasta"

		for align_type in $align_types; do

			#This is slightly redundant as we're doing it for each loop.
			ValidateVariable align_type VALID_ALIGN_TYPES
			seqs_var="${align_type}SEQS" 
			CheckVariables $seqs_var
			
			if [ $align_type = GENOMIC ]; then
				logic_name="Submit_${format}_ProbeAlign"
				
			elif [ $align_type = TRANSCRIPT ]; then
				logic_name="Submit_${format}_ProbeTranscriptAlign"
			fi
		
					
			#We could tune the batch size here
			#But we're already tuning the chunks to ~30mins
			#So we know that this should complete within about 3 hours
			#as we have a batch size of 5.
			#This doesn't always seem to work
			#As I have seen 126 running jobs from 500 chunks when batch
			#size is 5.  This should have just been 100, so batch size has been reduced
			#somehow

			#Execute perl $PIPELINE_SCRIPTS/rulemanager.pl $PDB_SCRIPT_ARGS -analysis $logic_name
			#Do we need input_id_type...will this not be picked up?
			#-input_id_type PROBE_CHUNK

         	#append to array list for species???
        	#echo $ARRAY >> ${WORK_DIR}/array.list
		done
	done

	echo "Now do ProbeAlignReport when finished"
}



################################################################################
# Func      : AlignWait
# Desc      : Accumulator step to wait for all ProbeAlign jobs
# Arg[1..n] : None
# Return    : None 
# Exception : None
################################################################################


AlignWait(){
	echo ":: AlignWait $*"
	#add opts here for -force_accumulators and -once?

	cmdline="$EFG_PERL $PIPELINE_SCRIPTS/rulemanager.pl $PDB_SCRIPT_ARGS -force_accumulators -analysis Align_Wait"

	#echo ":: ImportWait $cmdline"

	Execute $cmdline
	
	#echo "Can we cat the AFFY and AFFY_ST arrays_nr files?"
	#As these should use the same rules?  This would mean running with one logic name for both
	#Just keep separate for now

	#Now cat the array names files output by ImportArrays
	#This is to enable any ProbeFeature/IDXref/UnmappedObject deletes

	echo ": Finished waiting...ready for RunTranscriptXrefs?"
}


ProbeAlignReport(){
	echo ":: ProbeAlignReport $*"
			
	formats=
	custom_names=
	number=

	usage='usage: ProbeAlignReport [ -f "array_format"]* [ -a "array_name" ]* [ -n number of most mapped trancripts to list defualt 20 ][-h help]'

	OPTIND=1

	while getopts ":t:f:a:h" opt; do
		case $opt in 
	        f  ) formats="$OPTARG $formats";;
		    a  ) custom_names="$OPTARG $custom_names";;
		    n  ) number=$OPTARG;;
		    \? ) echo $usage; exit 1;;
			h  ) echo $usage; return 0;;
		esac 
	done		

	SetArrayFormats $formats
	number=${number:=20}


 	for format in $ARRAY_FORMATS; do
		SetArrayNames $format
		
		#This sets custom_names to $ARRAY_NAMES also!
		#names=${custom_names:=$ARRAY_NAMES}
		if [ $custom_names ]; then
			names=$custom_names
		else
			names=$ARRAY_NAMES
		fi

	
		#Validate custom names(sub this)?
		valid_names=
		
		for name in $names; do
			
			valid=$(echo $ARRAY_NAMES | grep -E "^$name | $name | $name$|^$name$")
			
			if [[ $valid ]]; then
				valid_names="'$name' $valid_names"
			fi
		done
		
		valid_names=$(echo $valid_names | sed 's/ /, /g')
				
		echo ""

		if [[ ! $valid_names ]]; then
			echo ": Skipping $format ProbeAlignReport, no valid custom names found: $names($ARRAY_NAMES)"
		else
			report_file=${WORK_DIR}/ProbeAlign.$format.report
			echo ": Generating $format($valid_names) ProbeAlignReport: $report_file"
			BackUpFile $report_file

	
			
			#Total Probes Mapped and Features
			
			#These should only be more than 100 if we get transcript mappings
			#Up to a maximum of 200
			

			query="SELECT 'Mapped Probes' as ''; "
			#Need to -skip-column-names (-N) here
			#This breaks the other result table formatting
			query="$query SELECT rpad('Array Name', 20, ' '), rpad('Align Analysis', 40, ' '), rpad('Total Probes Mapped', 19, ' '), rpad('Total Probe Features', 20, ' '); "
			logic_name=

			#Changed this to a logic_name loop
			#as we could accurately count the nulls

			for align_type in $ALIGN_TYPES; do

				if [ $align_type = TRANSCRIPT ]; then
					logic_name="${format}_ProbeTranscriptAlign"
				else
					logic_name="${format}_ProbeAlign"
				fi

				query="$query SELECT rpad(a.name, 20, ' '), rpad('$logic_name', 40, ' '), rpad(concat(count(distinct pf.probe_id), '/', count(distinct p.probe_id)), 19, ' '), rpad(count(pf.probe_feature_id), 20, ' ') from array a, array_chip ac, probe p left join (probe_feature pf, analysis an) on (p.probe_id=pf.probe_id and pf.analysis_id=an.analysis_id and an.logic_name='$logic_name') where a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id and a.name in($valid_names) group by a.array_id;"

			done


		
			
			#Probe frequencies

			#We can't rpad the counts here as this will make the sort lexical!
			query="$query select ' '; SELECT 'Probe Frequencies' as ''; select rpad('Feature Count', 13, ' '), rpad('Probe Count', 11, ' ');"
			query="$query select rpad(t.cnt, 13, ' '), rpad(count(t.cnt), 11, ' ') from (select count(distinct pf.probe_feature_id) as cnt from probe_feature pf, probe p, array_chip ac, array a WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND pf.probe_id=p.probe_id AND a.name in ($valid_names) GROUP by pf.probe_id) as t GROUP by t.cnt order by t.cnt;"
	
		    #Top 20 most abundant probes
			query="$query select ' '; select 'Top $number Most Mapped Probes' as '';"

			#Don't order by rpad'd FeatureCount as this makes the sort lexical
			#would have to do subselect otherwise(as above)
			#This is assuming all probes are different within a probeset
			#For those which are NOT tsk tsk
			#We need to group
			
			query="$query select rpad('FeatureCount', 15, ' '), 'ProbeID';"
			query="${query} select rpad(count(distinct pf.probe_feature_id),15, ' ') as FeatureCount, pf.probe_id from probe_feature pf, probe p, array_chip ac, array a WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND pf.probe_id=p.probe_id AND a.name in ($valid_names) GROUP by pf.probe_id order by count(distinct pf.probe_feature_id) desc limit $number;"
			
			
			#Total unmapped probes
			query="$query select ' '; select rpad('Total Unmapped Probes', 21,' ');"
			query="$query select rpad(count(distinct uo.ensembl_id), 21,' '), ur.summary_description from  probe p, array_chip ac, array a, unmapped_object uo, unmapped_reason ur WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND uo.ensembl_id=p.probe_id AND a.name in ($valid_names) and uo.unmapped_reason_id=ur.unmapped_reason_id and (ur.summary_description='Unmapped probe' or ur.summary_description='Promiscuous probe') GROUP by ur.summary_description;"


			#Top 20 most unmapped probe
			query="$query select ' '; select 'Top $number Most Unmapped Probes' as '';"
			#Should group here but this seems to hide some records
			#This sorts lexically! on full description, no mapping count, would need order by replace(ur.full_description, '/%', '')?
			query="$query select distinct rpad('ProbeID', 10, ' '), 'Reason';"
			query="$query select distinct rpad(p.probe_id, 10, ' '), ur.full_description from probe p, array_chip ac, array a, unmapped_object uo, unmapped_reason ur WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND uo.ensembl_object_type='Probe' AND uo.ensembl_id=p.probe_id AND a.name in ($valid_names) and uo.unmapped_reason_id=ur.unmapped_reason_id and ur.summary_description='Promiscuous probe' order by ur.full_description desc limit $number;"


	        #Mapped Probes by seq_region name and array
			#This is too big, need to split into array/seq_region queries
			#Needs to be done in perl script
			#query="${query} select ''; select 'Mapped Probes by array.name seq_region.name' as '';"
			#query="${query} select a.name, count(distinct pf.probe_feature_id), sr.name from probe_feature pf, probe p, seq_region sr, array_chip ac, array a where a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id and p.probe_id=pf.probe_id and pf.seq_region_id=sr.seq_region_id and a.name in ($valid_names) group by a.name, sr.name';


			#-N supress column names
			echo $query | mysql -N $DB_MYSQL_ARGS > $report_file
			cat $report_file

		
			

		fi
	done
}
	

#Should need to routinely use this as we have some of this info logged
#May need to use this if the log has been lost or needs validating

ProbeXrefReport(){
	echo ":: ProbeXrefReport $*"
	

	echo "Not yet implemented. See logs from probe2transcript.pl for xref info"
	return
		
	formats=
	custom_names=
	number=

	usage='usage: ProbeAlignReport [ -f "array_format"]* [ -a "array_name" ]* [ -n number of most mapped trancripts to list defualt 20 ][-h help]'

	OPTIND=1

	while getopts ":t:f:a:h" opt; do
		case $opt in 
	        f  ) formats="$OPTARG $formats";;
		    a  ) custom_names="$OPTARG $custom_names";;
		    n  ) number=$OPTARG;;
		    \? ) echo $usage; exit 1;;
			h  ) echo $usage; return 0;;
		esac 
	done		

	SetArrayFormats $formats
	number=${number:=20}


 	for format in $ARRAY_FORMATS; do
		SetArrayNames $format
		
		#This sets custom_names to $ARRAY_NAMES also!
		#names=${custom_names:=$ARRAY_NAMES}
		if [ $custom_names ]; then
			names=$custom_names
		else
			names=$ARRAY_NAMES
		fi

	
		#Validate custom names(sub this)?
		valid_names=
		
		for name in $names; do
			
			valid=$(echo $ARRAY_NAMES | grep -E "^$name | $name | $name$|^$name$")
			
			if [[ $valid ]]; then
				valid_names="'$name' $valid_names"
			fi
		done
		
		valid_names=$(echo $valid_names | sed 's/ /, /g')
				
		if [[ ! $valid_names ]]; then
			echo ": Skipping $format TranscriptXrefReport, no valid custom names found: $names($ARRAY_NAMES)"
		else
			report_file=${WORK_DIR}/TranscriptXref.$format.report
			echo ": Generating $format($valid_names) TranscriptXrefReport: $report_file"
			BackUpFile $report_file

			#Total Xref'd probes
			#Total xrefs per array

	#select count(distinct ensembl_id), a.name from object_xref ox, xref x, probe p, array_chip ac, array a where a.class='AFFY_UTR' and a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id and p.probe_set_id=ox.ensembl_id and ox.ensembl_object_type='ProbeSet' group by a.name with rollup;
	#With rollup isn't appropriate here as probesets are resundant across arrays?

	#WARNING!! This seems to cause index failure on tmp tables, maybe because tmp table get's too big?


	#Do we need a count across seq_regions just for the whole class?
	#This would require join to transcript table, so not possible at present without doing some horrible sql on individual probe xrefs

#mysql> select count(*), ur.summary_description  from unmapped_object uo, unmapped_reason ur, probe p, array_chip ac, array a where a.class='CODELINK' and a.array_id=ac.array_chip_id and ac.array_chip_id=p.array_chip_id and p.probe_id = uo.ensembl_id and uo.ensembl_object_type='Probe' and uo.unmapped_reason_id=ur.unmapped_reason_id;



			#Top 5 most Xref'd transcripts
			#Top 5 most Xref'd probes(inc prom)
			#Top 5 most Xref'd probes(no prom)
			
			#Total Probes Mapped and Features
			#Is this correct?
			#Why are we getting >300 freqs here?
			#These should only be more than 100 if we get transcript mappings
			#Up to a maximum of 200
			query="SELECT 'Annotated Probes' as ''; SELECT a.name as 'Array Name', an.logic_name as 'Align Analysis', count(distinct pf.probe_id) as 'Total Probes Mapped', count(pf.probe_feature_id) as 'Total Probe Features' from analysis an , array a, array_chip ac, probe p, probe_feature pf where a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id and p.probe_id=pf.probe_id and pf.analysis_id=an.analysis_id and a.name in($valid_names) group by a.array_id, pf.analysis_id;"
	
			
			#Probe frequencies
			query="${query} SELECT 'Probe Frequencies' as ''; select t.cnt as 'Feature Count', count(t.cnt) as 'Probe Count' from (select count(pf.probe_id) as cnt from probe_feature pf, probe p, array_chip ac, array a WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND pf.probe_id=p.probe_id AND a.name in ($valid_names) GROUP by pf.probe_id) as t GROUP by t.cnt order by t.cnt;"
	
		    #Top 20 most abundant probes
			query="${query} select ''; select 'Top $number Most Mapped Probes' as '';"
			query="${query} select count(pf.probe_id) as FeatureCount, pf.probe_id as ProbeID from probe_feature pf, probe p, array_chip ac, array a WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND pf.probe_id=p.probe_id AND a.name in ($valid_names) GROUP by pf.probe_id, a.array_id order by FeatureCount desc limit $number;"
	
				#Total unmapped probes
			query="${query} select ''; select count(uo.ensembl_id) as 'Total Unmapped Probes', ur.summary_description from  probe p, array_chip ac, array a, unmapped_object uo, unmapped_reason ur WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND uo.ensembl_id=p.probe_id AND a.name in ($valid_names) and uo.unmapped_reason_id=ur.unmapped_reason_id and (ur.summary_description='Unmapped probe' or ur.summary_description='Promiscuous probe') GROUP by ur.summary_description;"


			#Top 20 most unmapped probe
			query="${query} select ''; select 'Top $number Most Unmapped Probes' as '';"
			#Should group here but this seems to hide some records
			query="${query} select p.probe_id as 'Probe ID', ur.full_description as 'Reason' from probe p, array_chip ac, array a, unmapped_object uo, unmapped_reason ur WHERE ac.array_id=a.array_id AND ac.array_chip_id=p.array_chip_id AND uo.ensembl_object_type='Probe' AND uo.ensembl_id=p.probe_id AND a.name in ($valid_names) and uo.unmapped_reason_id=ur.unmapped_reason_id and ur.summary_description='Promiscuous probe' order by ur.full_description desc limit $number;"
			echo $query | mysql $DB_MYSQL_ARGS > $report_file
		fi
	done
}
	

#Do UnmappedObjectReport after probe2transcript


##Recover AWOL jobs?
#Lists disinct states of AWOL jobs which have had there AWOL status deleted?
#select distinct status from job_status where job_id not in(select j.job_id from job_status js, job j where j.job_id=js.job_id and is_current='y');

#Not sure about this as there is no SUCCESFUL state?
#Are jobs simply removed from table when successful?


#Change this to RunProbe2Transcript
#Then have RunTranscriptXrefs as a wrapper which would also call TranscriptXrefsWait, tail -n 25 each log and GenerateCDF etc...

#Simply sum the files used and multiply by factor


SetMusage(){
	array_file=
	factor=
	type=
	base_mem=

	OPTIND=1
	
	usage='write usage!'

	while getopts ":a:f:t:b:h" opt; do
		case $opt in 
	        a  ) array_file=$OPTARG;;
		    f  ) factor=$OPTARG;;
		    t  ) type=$OPTARG;;
			b  ) base_mem=$OPTARG;;
		    \? ) echo $usage; exit 1;;
			h  ) echo $usage; return 0;;
		esac 
	done		

	#Do we need to check factor is num?
	CheckVariables array_file factor

	seq_file_type=
	total_size=0

	if [ $type ]; then
		seq_file_name="${type}SEQS"
		CheckVariables $seq_file_name
		seq_file=$(eval "echo \$$seq_file_name")
		seq_file_type=seq
		seq_size=
	fi


	for file_type in $seq_file_type array; do
		file=$(eval "echo \$${file_type}_file")
		SetFileSize $file		
		total_size=$(($file_size + $total_size))
	done


	#This should work but doesn't seem to!!

	#if [ $type ]; then
	#	seq_file_name="${type}SEQS"
	#	CheckVariables $seq_file_name
	#	file=$(eval "echo \$$seq_file_name")
	#	echo "seq_file is $seq_file"
	#	total_size=$(GetFileSize $file)
		#upweight the array_file by cubing?
		#Is this valid for all the memusage calcs?
		#No this brings values closer!??
		#How is this possible if AFFY_ST is larger, growth should be exponential???

		#total_size=$(echo "($total_size * $total_size) / $WEIGHT" | bc) 
		#total_size=$(echo "($total_size * $total_size * $total_size) / $WEIGHT" | bc) 
		#total_size=$(echo "($total_size * $total_size * $total_size * $total_size * $total_size) / $WEIGHT" | bc) 
	#fi

	#eval ${file_type}_size=$(GetFileSize $file)
	#size=$(GetFileSize $array_file)
	#total_size=$(($size + $total_size))
	
	
		
	musage=$(echo "($total_size * $factor)/1000" | bc)
    musage=$(echo $musage | sed 's/\..*//')

	#Add base usage
	if [ $base_mem ]; then
		musage=$(($musage + $base_mem))
	fi


}



#HUMAN_AFFY_ST_probe2transcript select[mem>23256] rusage[mem=23256]



#Simply make a working dir

MakeHugememWorkdir(){
	echo ": MakeHugememWorkdir $*"

	format=$1
	shift

 	hugemem_workdir="${HUGEMEM_HOME}/$DB_NAME"
	echo "making $hugemem_workdir"

	if [ ! -d $hugemem_workdir ]; then
		Execute mkdir -p $hugemem_workdir
	fi

	#Create links here to err log out files?
	#Can we do this before they have been created?
	#This is step specific so leave to callers

}


RunTranscriptXrefs(){
    #$1 is update flag, i.e. this will delete all xrefs and object_xrefs
    #Speed up update process
	echo ":: RunTranscriptXrefs $*"
	

	formats=
	delete_xrefs=
	array_names=
	skip_backup=
	local=
	no_execute=
	memory=
	debug=
	usage='usage: RunTranscriptXrefs [ -f "array_format" default=$ARRAY_FORMATS ]+ [ -d(elete xrefs) -s(kip backup) -h(elp) -l(ocal i.e. no bsub) -n(o execute i.e. just echo cmd) -m(emory in kb if array_nr fasta not available) -b (debug mode)]'
	
	OPTIND=1
	
	while getopts ":f:a:m:lnbsdh" opt; do
		case $opt in 
	        f  ) formats="$OPTARG $formats";;
		    a  ) array_names="$OPTARG $array_names";;
		    d  ) delete_xrefs="--delete";;
	        s  ) skip_backup=1;;
			l  ) local=1;;
			m  ) memory=$OPTARG;;
			n  ) no_execute=1;;
            b  ) debug='--debug';;
		    \? ) echo $usage; exit 1;;
			h  ) echo $usage; return 0;;
		esac 
	done		

	#Delete is not being picked up on first call, but after??
	#Problem was var name clash in GenerateArrayNames, now fixed but still risk with other funcs due to lack of scoping?

	if [ "$QUEUE_MANAGER" = Local ]; then		
		local=1;
	fi

	
	echo ": BSUB_PROBE2TRANSCRIPT is $BSUB_PROBE2TRANSCRIPT"
	echo ": PROBE2TRANSCRIPT_PARAMS are $PROBE2TRANSCRIPT_PARAMS"
	

	SetArrayFormats $formats
	pwd=$PWD
    ChangeMakeDir $WORK_DIR
	#Move to wrapper when we create one
	if [ ! $skip_backup ]; then
		BackUpTables -t xrefs -s orig
	fi
	
  
	#Need to optionally take array names, as we may not have done the Alignment with this DB
	#i.e. if we are just rerunning the Xrefs due to a new genebuild
    #if [[ ! -e ${WORK_DIR}/arrays.list ]]; then echo " no arrays.list file error here!!!"; return 101; fi

    #check meta_coord table
    echo "Checking meta_coord table"
	#We need to restrict this to the specific analyses
    coord_id=$(QueryVal OUT "select coord_system_id from probe_feature pf, seq_region sr where pf.seq_region_id=sr.seq_region_id limit 1");
    meta_ids=$(QueryVal OUT "select coord_system_id from meta_coord where table_name='probe_feature'")

    cnt=0

	if [ $coord_id ]; then

		for id in $meta_ids; do

			if [ $id -eq $coord_id ]; then
       		#echo "Found correct level entry for oligo_feature coord_system"
				cnt=`expr $cnt + 1`
			fi
			
		done
	fi

    if [ $cnt -eq 0 ]; then
		echo "Cannot find valid meta_coord entry for probe_feature, maybe you have nor probe_features or did you migrate the tables from a different DB?"
		exit 0;
       	#echo "Generating meta_coord entry for oligo_feature with coord_id $coord_id"
       	#echo "insert into meta_coord values ('probe_feature', $coord_id, 25)" | mysql $DB_MYSQL_ARGS
    fi

	#No healtcheck done anymore as we just delete anyway and fail if checks fail
	#??? The healthcheck was to prevent deleting xrefs if we have associated arrays?

	
	for format in $ARRAY_FORMATS; do
		#Check we have an array.names file
		#To see whether we have imported
		NAMES_FILE="$WORK_DIR/arrays.${format}.names"

		if [[ ! -f $NAMES_FILE ]]; then
			echo "NAMES_FILE not found: $NAMES_FILE"
			GenerateArrayNames -f $format
		fi

		echo "Running $format TranscriptXrefs"
		SetArrayNames $format


		if [ $array_names ]; then
			names=$array_names
		else
			names=$ARRAY_NAMES
		fi
		
		if [[ ! $names ]]; then
			echo "WARNING: No $format ARRAY_NAMES found in: $NAMES_FILE"
			echo "Maybe you want to SetArrayNames $format ...?"
			exit 1;
		else

			workdir


				#Validate names
			valid_names=

			for name in $names; do
				
				valid=$(echo $ARRAY_NAMES | grep -E "^$name | $name | $name$|^$name$")
				
				if [[ $valid ]]; then
					valid_names="$name $valid_names"
				fi
			done
			
			
			if [[ ! $valid_names ]]; then
				echo "WARNING: No valid custom ARRAY_NAMES: $names($ARRAY_NAMES)"
				echo "A typos in your custom array names or maybe?"
				echo "Or maybe you want to resitrcit the ARRAY_FORMATS by using: SetArrayFormats format?"
				exit 1;
			fi

			
			vendor=$(_getVendorByFormat $format)
				
			if [ $? -ne 0 ]; then
				echo -e $vendor
				exit 1;
			fi

				#Here we need options to pass params in case format isn't supported by probe2transcript #Also need to calculate memusage by counting nr_fasta #What if not present
					        #echo "bsub cmd is $BSUB_CMD"
		        #Check here to see if queue is hugemem and exit if EFG_DATA does not begin with $HOME?
		        #Or just warn anyway
			BSUB_CMD=

			if [ -n "$BSUB_PROBE2TRANSCRIPT" ]; then
				file_name="${format}_probe2transcript"
				job_name="${DB_NAME}_${file_name}"
					#set musage
		
				if [ $memory ]; then
					musage_k=$memory
					musage=$(($memory/1000))
				else
					SetMusage -a "$WORK_DIR/arrays_nr.${format}.fasta" -f $XREF_MEM_FACTOR -t TRANSCRIPT -b $XREF_MEM_BASE
					musage_k=$(echo "$musage * 1000" | bc)
				fi
				
				work_dir=

				if [ $musage -gt $MAX_HUGE_MEM ]; then
					echo "Expected memory usage exceeds that available on $HUGEMEM_QUEUE queue: $musage. Modify rsuage generation?"
					exit 1;
				elif [ $musage -gt $MAX_NORMAL_MEM ]; then
					queue=$HUGEMEM_QUEUE


						#Is this a pipeline method?
						#Do we have to copy anything back afterwards?
					MakeHugememWorkdir $format
					work_dir=$hugemem_workdir

						#Now create links to files which will be created
					link_names="${file_name}.err ${file_name}.out ${DB_NAME}_${format}_probe2transcript.log ${DB_NAME}_${format}_probe2transcript.out"
						
					for link_name in $link_names; do
							
						if [ ! -L $link_name ]; then
							Execute ln -s "${work_dir}/${link_name}" $link_name
						fi
						
					done
				else
					queue=$VLONG_QUEUE
					work_dir=$WORK_DIR
				fi

					#We cannot let lsf checkpoint as this may be in the middle of a query
					#So we check point in the script pointing to the file name lsf would use
					#And then set the check point time greater than the time limit of the queue
					#But does LSF checkpoint automatically before exiting?  Which would overwrite the checkpoint file!

					#BSUB_CMD="-k \"${WORK_DIR}/checkpoints 1 method='blcr'\" -q $queue -R\"select[mem>$musage] rusage[mem=$musage]\" -M $musage_k -e ${work_dir}/${file_name}.err -o ${work_dir}/${file_name}.out cr_run " # -J $job_name"removed as this is now done in submitJob

				BSUB_CMD="-q $queue -R\"select[mem>$musage] rusage[mem=$musage]\" -M $musage_k -e ${work_dir}/${file_name}.err -o ${work_dir}/${file_name}.out"
			fi
			
			transpass=

			if [ $DNADB_PASS ]; then
				transpass="--transcript_pass $DNADB_PASS"
			fi

			script_cmd="$EFG_PERL ${EFG_SRC}/scripts/array_mapping/probe2transcript.pl --species $SPECIES --transcript_dbname $DNADB_NAME --transcript_host $DNADB_HOST --transcript_port $DNADB_PORT --transcript_user $DNADB_USER $transpass --xref_host $DB_HOST --xref_dbname $DB_NAME --xref_user $DB_USER --xref_pass $DB_PASS $PROBE2TRANSCRIPT_PARAMS $delete_xrefs $MULTI_SPECIES -vendor $vendor -format $format -arrays $valid_names -import_edb $debug"
			
			
			if [ -n "$BSUB_PROBE2TRANSCRIPT" ] && [ ! $local ]; then
					#bjobs -w -J $job_name
				
				if [ ! $no_execute ]; then
					submitJob $job_name "$BSUB_CMD" "$script_cmd"
				else
					echo "Skipping job"
					echo  $job_name $BSUB_CMD $script_cmd
				fi		
				
			else
			
				if [ ! $no_execute ]; then
					echo $script_cmd
					Execute $script_cmd
					echo "ps -ef | grep 'probe2transcript' to find out if it is still running"
					echo "Alternatively tail -f $WORK_DIR/${DB_NAME}_${format}_probe2transcript.log"
				else
					echo "Skipping job"
				fi
				
				echo "$script_cmd"	
				
			fi
			
		fi
		
	done
	
	cd $pwd
#echo "Remember to run TranscriptXrefsWait, GenerateCDF"
}

	

#This is a work around method for formats
#Which do not following the naming convention
#format name = VENDOR_FORMAT e.g. AFFY_UTR
#Formats which break this convention include
#LEIDEN, vendor = SPAINK_LAB_LEIDEN
#This is basically to avoid having to set vendor values 
#for each array format

_getVendorByFormat(){
	format=$1

	if [ ! $format ]; then
		echo "You must supply a format argument to _getVendorByFormat"
		return 1
	fi

	vendors=($(QueryVal OUT "select distinct vendor from array where class='$format'"))

	if [[ ${#vendors[*]} -gt 1 ]]; then
		echo -e "Could not find unique vendor for $format format:\t$vendors"
		echo -e "Please specify valid vendor using -v option"
		return 1
	else
		vendor=$vendors

		if [ ! $vendor ]; then
			echo "Could not find $format in the DB, maybe you have not yet imported this array format?"
			return 1
		fi
	fi

	echo $vendor
}


### Miscellaneous methods

#Move these to base env Get***by_DB?
#And then have these call the efg.env method
#So we don't have to sourve the arrays.env so use them?
#Add UnmappedObject, Probe and ProbeFeature methods
#Or write a quick tools script and have these methods call the tools script
#Then we can have this func in the API too
#Write wrapper run script which has appropriate help and uses default env vars 
#Then these methods can just call the run script

GetTranscriptXrefs(){
	echo ": GetTranscriptXrefs $*"
	
	sids=
	usage='usage: GetTranscriptXrefs [ -s TRANSCRIPT_STABLE_ID ]+ [ -h(elp) ]'
	
	OPTIND=1
	
	while getopts ":s:h" opt; do
		case $opt in 
	        s  ) sids="$OPTARG $sids";;
		    \? ) echo $usage; return 1;;
			h  ) echo $usage; return 0;;
		esac 
	done		
	
	CheckVariablesOrUsage "$usage" sids				

	sids=$(echo $sids | sed "s/ /', '/")	


	#2nd arg as TARGET|OUT
	query="select ps.name, x.dbprimary_acc, ac.name, edb.db_release from probe_set ps, array_chip ac, probe p, xref x, object_xref ox, external_db edb where edb.external_db_id=x.external_db_id and x.dbprimary_acc in('$sids') and ps.probe_set_id=ox.ensembl_id and ox.ensembl_object_type='ProbeSet' and ox.xref_id=x.xref_id and ps.probe_set_id=p.probe_set_id and p.array_chip_id=ac.array_chip_id group by ps.name order by x.dbprimary_acc, ac.name, ps.name, edb.db_release"


	echo $query | mysql $MYSQL_ARGS

}


GetProbesetXrefs(){
	echo ": GetProbesetXrefs $*"
	
	names=
	arrays=
	usage='usage: GetProbesetXrefs [ -n name ]+ [ -a array_name ]* [ -h(elp) ]'
	
	OPTIND=1
	
	while getopts ":n:a:h" opt; do
		case $opt in 
	        n  ) names="$OPTARG $names";;
	        a  ) arrays="$OPTARGS $arrays";;  
		    \? ) echo $usage; return 1;;
			h  ) echo $usage; return 0;;
		esac 
	done		
	
	CheckVariablesOrUsage "$usage" names			

	names=$(echo $names | sed "s/ /', '/")	
	arrays=$(echo $arrays | sed "s/ /', '/")

	array_clause

	if [ $arrays ]; then
		array_clause=" and ac.name in('$arrays') "
	fi
	



	#2nd arg as TARGET|OUT
	query="select ps.name, x.dbprimary_acc, ac.name, edb.db_release from probe_set ps, array_chip ac, probe p, xref x, object_xref ox, external_db edb where edb.external_db_id=x.external_db_id and ps.name in('$names') and ps.probe_set_id=ox.ensembl_id and ox.ensembl_object_type='ProbeSet' and ox.xref_id=x.xref_id and ps.probe_set_id=p.probe_set_id and p.array_chip_id=ac.array_chip_id $array_clause group by ps.name order by x.dbprimary_acc, ac.name, ps.name, edb.db_release"


	echo $query | mysql $MYSQL_ARGS

}



#Is this not in funcs?

ContinueOverride(){
    CheckVariables 1

    if [ ! $2 ]
    then
		AskQuestion "$1"

		if [[ "$REPLY" != [yY]* ]]
		then
			echo "Exiting"
			exit
		fi
    else
		echo "Auto Continue"
    fi

}

CompareStagingMirrorDBs(){
	echo "This method needs updating to use both staging servers"
	return

    SVER=$1
  
    CheckVariables SVER

    echo "Comparing assembly/genebuild infro between staging and live servers"

    MVER=`expr $SVER - 1`


    STAGING="-hens-staging -uensro"
    ENSDB="-hensembldb.ensembl.org -uanonymous"
    CNT=0

    MYSQL_ARGS=$STAGING
    SDBS=$(QueryVal show databases like \"%core_${SVER}_%\")

    for dbname in $SDBS
    do
	CNT=`expr $CNT + 1`

	if [ $CNT -gt 2 ]
	then
      	    MYSQL_ARGS="$STAGING $dbname"
	    GBUILD[$CNT]=$(QueryVal select meta_value from meta where meta_key=\"genebuild.version\")
	    ASSEMBLY[$CNT]=$(QueryVal select meta_value from meta where meta_key=\"assembly.name\")
	fi
    done

    CNT=0
    MYSQL_ARGS=$ENSDB


    OFILE="${AFFY_HOME}/assembly_build_comparison.${SVER}_${MVER}"

    if [ -f $OFILE ]
    then
	rm -f $OFILE
    fi


    for dbname in $SDBS
    do
	CNT=`expr $CNT + 1`

	if [ $CNT -gt 2 ]
	then
	    edbname=$(echo $dbname | sed s"/_core_${SVER}_[0-9]*[a-z]*/_core_${MVER}/")   
      	    MYSQL_ARGS="$ENSDB"
	    edbname=$(QueryVal show databases like \"${edbname}%\");

	    if [[ $edbname ]]
	    then
		OUTLINE=""
		edbname=$(echo $edbname | sed 's/.*)//')
		MYSQL_ARGS="$ENSDB $edbname"
		EGBUILD=$(QueryVal select meta_value from meta where meta_key=\"genebuild.version\")
		EASSEMBLY=$(QueryVal select meta_value from meta where meta_key=\"assembly.name\")
	      
		#compare dbnames here to get genebuild difference
		EGB=$(echo $edbname | sed s'/.*_core_[0-9]*_//')
		SGB=$(echo $dbname | sed s'/.*_core_[0-9]*_//')

		if [[ ${ASSEMBLY[$CNT]} != $EASSEMBLY ]]
 		then
		    OUTLINE="Mapping ({$EASSEMBLY} > ${ASSEMBLY[$CNT]})  &  "
		fi


		if [[ ${GBUILD[$CNT]} != $EGBUILD ]]
 		then
       		    OUTLINE="$OUTLINE    Xrefing (${EGBUILD} > ${GBUILD[$CNT]})"
		elif [[ $EGB != $SGB ]]
		then
		    OUTLINE="$OUTLINE    Xrefing (Patched set $EGB > $SGB contact genebuilder?)"
		elif [[ $OUTFILE ]]
		then
		    OUTLINE="$OUTLINE  WARNING: SAME GENEBUILD FOR NEW ASSEMBLY"
		fi

	  
		if [[ $OUTLINE ]]
		then
		    OUTLINE="$dbname requires:    ${OUTLINE}"
		    echo $OUTLINE
		    echo $OUTLINE >> $OFILE
		fi


	    else
		echo "New DB $dbname!!!"
	    fi
	fi
    done

    #get all names from staging matching core_SVER
    #for each add element to assembly and genebuild arrays
    
    #for each one, sed for core_MVER
    #then test assembly and genebuild vs array values
    

}


CountProbesetsPerArray(){
	arrays=$*

	#Validate array exists here by grabbing array_chip_id
	
	in_arrays=

	for array in $arrays; do
		in_arrays="'${array}', "
	done

	in_arrays=$(echo $in_arrays | sed 's/,$//')



	sql="select a.name, count(distinct probe_set_id) from probe p, array_chip ac, array a where a.name in($in_arrays) and a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id group by a.name"

	echo $sql

	mysqlefg -e "$sql"

}

#Largely redundant as ProbeAlign checks for IMPORTED status

ShowImportedArrays(){

	sql="select rpad(ac.name, 30, ' ') as 'ArrayChip\t\t', rpad(a.class, 20, ' ') as 'Class\t\t', rpad(s.name, 15, ' ')  as 'Imported Status' from array a, array_chip ac left join (select s1.table_id, sn1.name from status s1, status_name sn1 where s1.status_name_id=sn1.status_name_id and s1.table_name='array_chip' and sn1.name='IMPORTED') s on s.table_id=ac.array_chip_id where a.array_id=ac.array_id and a.class is not NULL order by a.class, a.name"

	echo $sql | mysqlefg

}
