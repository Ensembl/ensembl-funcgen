#!/usr/local/bin/bash

#Test for eFG env
if [  -z "$EFG_SRC" ] || [ ! -d $EFG_SRC ]; then
   echo ":: You have not yet initialised the eFG environment"
   return	
fi

#still have to source the funcs!
#As not exported in env?
. $funcs_file


#To do
#1 Integrate this into pipeline.env?
#  This requires set up of a instance file and input dir. Also is less visible? 
#  Need to fettle with Help to list funcs better?
#  As we have no structure to the functions, but at least we have a dir 
#  structure when we have separate scripts, add outdir ?
#2 remove hardcoded fastagrep path
#3 Done redirect output to standardised alignments output directory structure
#4 Done Take any input path, so we don't have to have standardised structure for fastq dirs as this often mirrors source repository structure
#5 Add -cell_type and -feature_type opts and cat to CellType_FeatureType(_.*).fastq, or just add set_name?
#6 DONE Finish files handling
#7 DONE Add fastq.list file to log the true uncat'd input fastq file(s)  
#8 Optionally wait for job success and remove cat'd file if present?
#9 Use $EFG_DATA for scratch_dir.  Currently set to group rather than personal scratch dir
#  Need to handle both $EFG_DATA and $EFG_GROUP_DATA?
#10 DONE Integrate Daniel's batch code
#11 Wrap bsub cmd to catch errors and exit properly
#12 Add no clean up flag to keep intermediate files? 

OPTIND=1    #Reset the getopts if we have used it previously
VALID_GENDERS='male female'
gender=
species=
index_home=/lustre/scratch103/ensembl/funcgen/bwa_indexes    #Set this in pipeline/efg.config?
assembly=GRCh37
mask='_unmasked'
align_type='samse'
resource='-R"select[mem>5000] rusage[mem=5000]" -M5000000'
dir=
clean=
experiment_name=
name=
outdir=
format=sam
merge_only=0
#Change this to $EFG_DATA

scratch_dir=/lustre/scratch103/ensembl/funcgen

usage="Usage:\t\trun_bwa.sh <options> [ file1 file2 ]
\n\nDescription:\tThis script will take either a fastq file or a directory as input, unzip cat/rename 
\n\t\tand rezip as necessary to generate a single fastq file. The relevant bwa indexes will be validated 
\n\t\tbefore submitting the bwa job to the farm.
\n\nExample:\trun_bwa.sh  -g male -a GRCh37 -s homo_sapiens -e Stamatoyannopoulos_DNase1_EpiRoadmap_GSE18927  -n IMR90_DNase1 SRX01242*/*fastq*
\n\nMandatory options:
\n\t-s(pecies e.g. homo_sapiens)
\n\t-g(ender e.g male|female)
\n\t-n(set name i.e. CellType_FeatureType(_.*), in line with expected input for peaks.env)
\n\nOther options:
\n\t-e(xperiment name, required if input files are split over different directories)
\n\t-o(utdir overrides default $scratch_dir/alignments/\$species/\$experiment_name)
\n\t-a(ssembly default is $assembly)
\n\t-p(aired default is single reads)
\n\t-i(ndex_home default=$index_home)
\n\t-r(esource default=$resource)
\n\t-c(lean away cat'd fastq file)
\n\t-m(erge only, only submit the merge job if it has previously failed due to the absence of the same header)
\n\t-b(am output default=sam)
\n\t-h(elp)"


#Add opts:
#ftype
#ctype
#experiment name

#build alignments dir based on $EFG_DATA/alignments/species/assembly/exp_name
#remove file option and just have as args?


while getopts ":g:s:e:n:d:o:i:r:a:mbcph" opt; do
	case $opt in 
	    g  ) gender=$OPTARG ;; 
	    s  ) species=$OPTARG ;;
	    e  ) experiment_name=$OPTARG ;;
	    n  ) name=$OPTARG ;;
	    d  ) dir=$OPTARG ;;
	    o  ) outdir=$OPTARG ;;
	    i  ) index_file=$OPTARG ;;
	    r  ) resource=$OPTARG ;;
	    a  ) assembly=$OPTARG ;;
	    p  ) align_type='sampe' ;;
	    c  ) clean=1 ;;
	    m  ) merge_only=1 ;;
	    b  ) format='bam' ;;
	    h  ) echo -e $usage; exit 0;;
	    \? ) echo -e $usage; exit 1;;
	esac 
done


#Get trailing file name arguments
idx=1

while [[ $idx -lt $OPTIND ]]; do 
	idx=$(($idx + 1))
	shift 
done

files=($*)

#Parameter checks
CheckVariablesOrUsage "$usage" gender species index_home name
ValidateVariableOrUsage "$usage" gender VALID_GENDERS

if [[ ! "$files" ]] && [[ ! $dir ]]; then
	echo -e "ERROR:\tYou must supply either a -d(ir) parameter or some file arguments\n"
	echo -e $usage
	exit
elif [[ "$files" ]] && [[ $dir ]]; then
	echo -e "ERROR:You have specified -d(irectory) and input files, please use one or the other\n"
	echo -e $usage
	exit
fi


file_prefix=
new_input=1
unzipped_files=
split_files=


if [[ $dir ]] && [[ ! -d $dir ]]; then
    echo "-d $dir is not a valid directory"
    exit
fi

if [ ! "$experiment_name" ]; then
    if [[ "$files" ]]; then
	echo -e "ERROR:\tYou must explicitly set an -e(experiment name) when specify multiple file paths\n"
	echo -e $usage
	exit
    fi
    
    if [ ! "$dir" ]; then
	#Get exp name from dir name
	experiment_name=$(echo $dir | sed -r 's/\/$//')
	experiment_name=$(echo $experiment_name | sed -r 's/.*\///')
    fi
    
fi

if [ ! "$outdir" ]; then
    outdir="${scratch_dir}/alignments/${species}/${assembly}/${experiment_name}"
fi

file_prefix="${outdir}/${name}."
split_file="${file_prefix}1.fastq"

#Test for pre-cat'd file
if [[ -f $split_file ]]; then

    if [ $clean ]; then
	echo -e "Removing previously cached fastq files:\n\t"
	ls ${file_prefix}[1-9]*.fastq
	rm -f ${file_prefix}[1-9]*.fastq
    else 
	new_input=0
	split_files=($(ls ${file_prefix}[1-9]*.fastq))
	tmp=$(echo "${split_files[*]}" | sed 's/ /\n\t/g')
	echo -e "WARNING:\tUsing previously cached fastq files:\n\t${tmp}\n"
	echo -e "Specify -c(lean) to override this behaviour"
	sleep 5
    fi
else
    #Remove just incase we have other cached files hanging around
    rm -f ${file_prefix}[1-9]*.fastq
fi


#Get files
if [[ $new_input = 1 ]] && [[ $dir ]]; then
    files=($(ls ${dir}/*fastq*))
    
    if [[ ! "$files" ]]; then
	echo -e "ERROR:\tNo fastq files found in input dir:\t$dir\n"
	echo -e $usage
	exit
    fi
fi

#Do we not have a function for this unzipping/cating?

if [[ $new_input = 1 ]]; then
	zipped_files=

	for f in ${files[*]}; do
	
		if [[ $f != *fastq* ]]; then
			echo -e "Ignoring non-fastq file:\t$f"
		else
			CheckFile $f
			echo -e "Using fastq file:\t$f"
			funzipped=$f
			
	 	        #Unzip those which are zipped
			#would be nice to list to files first before unzipping, so we know what we are running with immediately


			if [[ $f =  *.gz ]]; then
				zipped_files=($f ${zipped_files[*]})
				funzipped=$(echo $f | sed 's/\.gz//')
			fi 
			
			unzipped_files=(${unzipped_files[*]} $funzipped)
		fi

	done
	
	if [ "$zipped_files" ]; then
		echo -e "Unzipping files..."	
		Execute gunzip ${zipped_files[*]}
	fi

	
	if [ ! "$unzipped_files" ]; then
		echo -e "ERROR:\tCould not find any valid input fastq files from list:\n${files[*]}"
		exit
	fi


fi

#Cat if required and log inputs
alignment_log="${outdir}/${name}.alignment.log"
input_file="$file_prefix\$LSB_JOBINDEX.fastq"



if [[ $new_input = 1 ]]; then
    split=1
    
    #Run from single input file		
    #Run from source if it is < 8000000 line long
    #This will preserve space on device and prevent uneccessary cat | split
	
    if [[ ${#unzipped_files[*]} = 1 ]]; then
		    
       if [[ $(cat ${unzipped_files[0]} | wc -l) -lt 8000000 ]]; then
    	  #cat | wc to avoid filename in output
	  echo "Skipping file batch as single input file is < 8000000 lines long"
	  split=0
	  echo "Input fastq file for $experiment_name $name alignments:" > $alignment_log
	  echo  "\t${unzipped_files[0]}" >> $alignment_log 
	  split_files=(${unzipped_files[0]})
	  input_file=${unzipped_files[0]}
       fi
    fi


    if [[ $split = 1 ]]; then
	MakeDirs $outdir			
	echo -e "Batching files to:\t${file_prefix}[1-9]\*\.fastq"
	#This is not catching failure!
	Execute cat ${unzipped_files[*]} | split -d -a 4 -l 8000000 - $file_prefix
	
	#rename files to add suffix, remove leading 0s
	#and +1 to the number to work with LSB_JOBINDEX values
	
	ls $file_prefix[0-9][0-9][0-9][0-9] | while read f; do num=$(echo $f | sed -r "s/.*\.([0-9]+)$/\1/"); num=$(echo $num | sed -r "s/^[0]+([0-9][0-9]*)/\1/"); num=$(($num + 1)); mv $f "${file_prefix}${num}.fastq"; done
	
	split_files=($(ls $file_prefix[1-9]*.fastq))
	echo -e "Created ${#split_files[*]} batch files"
	echo "Input fastq files for $experiment_name $name alignments:" > $alignment_log

	for f in ${unzipped_files[*]}; do
	    echo "\t${f}" >> $alignment_log
	done
	
        #Could maybe add a validation step here to wc -l the batches versus the input
        #echo -e "Total reads in fastq files" >> $alignment_log
        #echo -e $(expr $(cat ${unzipped_files[*]} | wc -l) / 4)  >> $alignment_log
	
	echo "Gzipping all source fastq files"
	Execute gzip ${unzipped_files[*]}
    fi
fi


#Should always have split_files by now
#but let's test for sanity
if [[ ! "$split_files" ]]; then
	echo -e "ERROR:\tCould not identify split_files. Something is wrong with the script!"
	exit
fi


#Could also log other run params here?
lc_species=$(echo $species | tr [A-Z] [a-z])
uc_species=$(echo $species | tr [a-z] [A-Z])
index_name="${lc_species}_${gender}_${assembly}${mask}.fasta"
fasta_file=${index_home}/${uc_species}/${index_name}

#Check index/fasta files
for ext in amb ann bwt pac rbwt rpac rsa sa; do
	error=$(CheckFile "${index_home}/${uc_species}/${index_name}.${ext}")

	if [ $? -ne 0 ]; then
		echo -e $error
		
		if [[ ! -f $fasta_file ]]; then
			echo -e "\nYou need to create the fasta file to index:\t$fasta_file"
			#Add how to dump here?
			echo -e "\ne.g. bsub -J ${lc_species}_top_level_unmasked -e ${lc_species}_top_level_unmasked.err -o ${lc_species}_top_level_unmasked.out -R\"select[mem>3500] rusage[mem=3500]\" -M3500000 $EFG_PERL $SRC/ensembl-analysis/scripts/sequence_dump.pl -dbhost ens-livemirror -dbuser ensro -dbname ${lc_species}_core_SCHEMA_BUILD -species $lc_species -coord_system_name toplevel -filename $lc_species_male_ASSEMBLY_unmasked.fasta"
			echo -e "Then create the female file as follows:\n\tcat ${lc_species}_male_ASSEMBLY_unmasked.fasta | /nfs/users/nfs_d/dkeefe/bin/fastagrep -t -v -P -p 'chromosome:ASSEMBLY:Y' > ${lc_species}_female_ASSEMBLY_unmasked.fasta"
		fi

		echo -e "\nYou need to generate the bwa indexes using:\n\tbsub -J ${index_name}_indexes -o ${index_name}_bwa_indexes.out  -e ${index_name}_bwa_indexes.err $resource bwa index -a bwtsw $fasta_file"
		
		exit
	fi
done



### Submit jobs
run_txt="\nRunning bwa with following options:\n\tIndex name\t= $index_name\n\tAlignment type\t= $align_type\n\tOutput dir\t= $outdir\n\tOutput format\t= $format\n"
echo -e $run_txt
echo -e $run_txt >> $alignment_log

### Run bwa and remove intermediate files...
#Use submitJob to avoid truncation of bsub cmd
align_job_name="bwa_${align_type}_${experiment_name}_${name}";
clean_cmd=

if [[ $merge_only = 1 ]]; then
	echo -e "Skipping alignment job:\t$align_job_name"
else

    #Is the index job still running?
	checkJob ${index_name}_indexes exit_if_running


	if [[ $new_data = 1 ]]; then
		echo -e "ERROR:\tCannot skip alignment job as there is new data, please omit -m(erge only)"
		exit;
	fi


	#we have problem with $LSB var being interpolated before they are set
	#Solution is to interpolate in the environment?
	#Do this passing '$job_cmd'
	#Will this still not interpolate when submitting job?

	bsub_cmd="-q long $resource -o ${outdir}/${align_job_name}.%J.%I.out -e ${outdir}/${align_job_name}.%J.%I.err"
	align_job_cmd="bwa aln $fasta_file  $input_file | "

	align_job_cmd="$align_job_cmd bwa $align_type $fasta_file - $input_file | "

#This bam sort is a little redundant at the moment as
#pipeline imports only take sam at present and always sorts before import.
#However, needed for merging and bam sort should be faster
#When we implement bam parsers, we can optionally remove the sam conversion
#and also add a sorted flag to the imports

	align_job_cmd="$align_job_cmd samtools view -uS - | "

	align_job_cmd="$align_job_cmd samtools sort - $file_prefix\$LSB_JOBID.\$LSB_JOBINDEX.${align_type}.sorted" 
	
	echo -e "\n"

	submitJob "$align_job_name[1-${#split_files[*]}]" "$bsub_cmd" '$align_job_cmd' eval
fi

#Now add merge and clean up job dependent on completion of first job
sam_header="${index_home}/${uc_species}/${lc_species}_${gender}_${assembly}${mask}.header.sam"

if [[ ! -f $sam_header ]]; then
	echo -e "ERROR:\tCould not find sam header to mere with:\t$sam_header"
	exit
fi

merge_cmd="samtools merge -h $sam_header - ${file_prefix}[0-9]*.[1-9]*.${align_type}.sorted.bam | "

merge_job_name="merge_${align_job_name}"

bsub_cmd=" -o ${outdir}/${merge_job_name}.%J.out -e ${outdir}/${merge_job_name}.%J.err"

if [[ $merge_only != 1 ]]; then
	bsub_cmd="$bsub_cmd -w 'done(${align_job_name}[1-${#split_files[*]}])' "
fi

# Consider removing duplicates...
#merge_cmd="$merge_cmd samtools rmdup -s - | " 

merge_cmd="$merge_cmd samtools view -h - | gzip -c > ${file_prefix}${align_type}.sam.gz"

clean_cmd=


echo "format is $format"

clean_cmd="rm -f ${file_prefix}[0-9]*.[1-9]*.${align_type}.sorted.bam"

# ouput some statistics
#possibly need to add the sam index (.fai file equivalent to header...)
#sam_index="${index_home}/${uc_species}/${lc_species}_${gender}_${assembly}${mask}.fa.fai"
log_cmd="echo \"Alignment QC - total reads as input: \" >> ${alignment_log}"
log_cmd="${log_cmd}; samtools view -uS ${file_prefix}${align_type}.sam.gz | samtools flagstat - | head -n 1 >> ${alignment_log}"
log_cmd="${log_cmd}; echo \"Alignment QC - mapped reads: \" >> ${alignment_log} "
log_cmd="${log_cmd}; samtools view -uS -F 4 ${file_prefix}${align_type}.sam.gz | samtools flagstat - | head -n 1 >> ${alignment_log}"
log_cmd="${log_cmd}; echo \"Alignment QC - reliably aligned reads (mapping quality >= 1): \" >> ${alignment_log}"
log_cmd="${log_cmd}; samtools view -uS -F 4 -q 1 ${file_prefix}${align_type}.sam.gz | samtools flagstat - | head -n 1 >> ${alignment_log}"
# Maybe do some percentages?

job_cmd="$merge_cmd; $clean_cmd; $log_cmd;"
#job_cmd="$merge_cmd; $clean_cmd;"

echo -e "\n"

submitJob "$merge_job_name" "$bsub_cmd" "$job_cmd"

#Could really do with submitting this to farm if we are waiting
#for cat gzip split activity
#Wait for job success here and remove cache file if present?
#Could also filter and convert sam2bed

echo -e "\nNOTE:\tNeed to remove the cat'd fastq file after the job has successfully completed"
