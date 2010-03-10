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

OPTIND=1    #Reset the getopts if we have used it previously
VALID_GENDERS='male female'
gender=
species=
file=
index_home=/lustre/scratch103/ensembl/funcgen/bwa_indexes    #Set this in pipeline/efg.config?
assembly=GRCh37
mask='_unmasked'
align_type='samse'
resource='-R"select[mem>5000] rusage[mem=5000]" -M5000000'
dir=
clean=
experiment_name=
outdir=
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
\n\t-h(elp)"


#Add opts:
#ftype
#ctype
#experiment name

#build alignments dir based on $EFG_DATA/alignments/species/assembly/exp_name
#remove file option and just have as args?


while getopts ":g:s:e:n:d:o:i:r:a:cph" opt; do
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






if [ "$dir" ]; then
	
	if [ ! -d $dir ]; then
		echo "-d $dir is not a valid directory"
		exit
	fi

    #Set outdir and experiment_name if not defined
	if [ ! "$experiment_name" ]; then
	    #Get exp name from dir name
		experiment_name=$(echo $dir | sed -r 's/\/$//')
		experiment_name=$(echo $experiment_name | sed -r 's/.*\///')
	fi

	if [ ! "$outdir" ]; then
		outdir="${scratch_dir}/alignments/${species}/${assembly}/${experiment_name}"
	fi

	
	file="${outdir}/${name}.fastq"

	#Test for pre-cat'd file
	if [[ -f $file ]]; then

		if [ $clean ]; then
			echo -e "Removing previously fastq file:\t${file}"
			rm -f $file
		else
			echo -e "WARNING:\tUsing previously cat'd fastq file:\t$file"
			echo -e "Specify -c(lean) to override this behaviour"
			files=($file)
			sleep 10
		fi
	fi


	#Get files
	if [[ ! -f $file ]]; then
		files=($(ls ${dir}/*fastq*))

		if [[ ! "$files" ]]; then
			echo -e "ERROR:\tNo fastq files found in input dir:\t$dir\n"
			echo -e $usage
			exit
		fi
	fi

else  # files

	#Test exp_name and set outdir and file
	if [ ! "$experiment_name" ]; then
		echo -e "ERROR:\tYou must explcitly set an -e(experiment name) when specify multiple file paths\n"
		echo -e $usage
		exit
	fi

	if [ ! "$outdir" ]; then
		outdir="${scratch_dir}/alignments/${species}/${assembly}/${experiment_name}"
	fi

	
	file="${outdir}/${name}.fastq"

	#Test for pre-cat'd file
	if [[ -f $file ]]; then

	 	if [ $clean ]; then
		 	echo -e "Removing previously cat'd fastq file:\t${file}"
			rm -f $file
		else 
 			echo -e "WARNING:\tUsing previously cat'd fastq file:\t$file"
			echo -e "Specify -c(lean) to override this behaviour"
	 		files=($file)
			sleep 10
	 	fi
 	fi
fi


#Do we not have a function for this unzipping/cating?
tmp_files=

for f in ${files[*]}; do
	
	if [[ $f != *fastq* ]]; then
		echo -e "Ignoring non-fastq file:\t$f"
	else
		CheckFile $f
		echo -e "Using fastq file:\t$f"
		
		funzipped=$f
		
	 	    #Unzip those which are zipped
		if [[ $f =  *.gz ]]; then
			echo -e "Unzipping..."

			exit

			Execute gunzip $f
			funzipped=$(echo $f | sed 's/\.gz//')
		fi 
		
		tmp_files=(${tmp_files[*]} $funzipped)
	fi
done

if [ ! "$tmp_files" ]; then
	echo -e "Could not find any valid input fastq files from list:\n${files[*]}"
	exit
fi

#Cat if required and log inputs
alignment_log="${outdir}/${name}.alignment.log"

#file up until now is the expected cat'd/cache fastq file 

if [[ ${#tmp_files[*]} -gt 1 ]]; then
	MakeDirs $outdir			
	echo -e "Cat'ing files to:\t$file"
	#This is not catching failure!
	Execute cat ${tmp_files[*]} > $file

	echo 'Input fastq files for $experiment_name $name alignments:' > $alignment_log
	
	for f in ${tmp_files[*]}; do
		echo $f >> $alignment_log
	done

	echo "Gzipping all source fastq files"
	Execute gzip ${tmp_files[*]}
else
	#Run from single input file

	#This may not be the full path? So match may fail
	#need to ls full patch for comparison

	if [[ ${tmp_files[0]} != $file ]]; then
		#Only log input when we know
		#we are not running with the cache file
		#otherwise preserve original inputs
		echo 'Input fastq file for $experiment_name $name alignments:' > $alignment_log
		echo  ${tmp_files[0]} >> $alignment_log
	fi

	file=${tmp_files[0]}
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


run_txt="\nRunning bwa with following options:\n\tIndex name\t= $index_name\n\tInput file\t= $file\n\tAlignment type\t= $align_type\n\tOutput dir\t= $outdir\n"
echo -e $run_txt
echo -e $run_txt >> $alignment_log

#Is the index job still running?
checkJob ${index_name}_indexes exit_if_running

#Use submitJob to avoid truncation of bsub cmd
#Could tidy up names a little here, no need 
set_name="${experiment_name}_${name}"
job_name="bwa_${align_type}_${set_name}";
bsub_cmd="-q long $resource -o ${outdir}/${job_name}.out -e ${outdir}/${job_name}.err"
job_cmd="'bwa aln $fasta_file  $file > ${outdir}/${set_name}.${align_type}.sai; bwa $align_type $fasta_file ${outdir}/${set_name}.${align_type}.sai $file > ${outdir}/${set_name}.${align_type}.sam'"

#echo $bsub_cmd $job_cmd

echo "\n"

submitJob "$job_name" "$bsub_cmd" "$job_cmd" 

#Could really do with submitting this to farm if we are waiting
#Wait for job success here and remove cache file if present?

echo -e "\nNOTE:\tNeed to remove the cat'd fastq file after the job has successfully completed"
