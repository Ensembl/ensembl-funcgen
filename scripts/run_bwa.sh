#!/usr/local/bin/bash



#Test for eFG env

if [  -z "$EFG_SRC" ] || [ ! -d $EFG_SRC ]; then
   echo ":: You have not yet initialised the eFG environment"
   return	
fi


#still have to source the funcs!
#As not exported in env
. $funcs_file


#To do
#Integrate this into pipeline.env?
#This requires set up of a instance file and input dir
#Also is less visible? Need to fettle with Help to list funcs better?
#As we have no structure to the functions, but at least we have a dir structure when we have separate scripts
#add code to generate indexes if absent?


#This makes sure we reset the getopts ind if we have used it previously
OPTIND=1
VALID_GENDERS='male female'
gender=
species=
file=
#Set this in pipeline/efg.config?
index_home=/lustre/scratch103/ensembl/funcgen/bwa_indexes
assembly=GRCh37
mask='_unmasked'
align_type='samse'
resource='-R"select[mem>5000] rusage[mem=5000]" -M5000000'


usage="usage:\trun_bwa.sh  -s(pecies e.g. homo_sapiens)  -g(ender e.g male|female)  -f(ile e.g. your_file.fastq)  [ -p(aired default is single reads)  -i(ndex_home default=$index_home)  -r(esource default=$resource)  -a(ssembly default=$assembly)  -h(elp) ] " #-m(asked default is unmasked)
#description='More wordy description here'

#To do
#add outdir 
#add inputdir option which will cat files? Or just use GetOptArgs on individual files?

while getopts ":g:s:f:i:r:a:ph" opt; do
	case $opt in 
		g  ) gender=$OPTARG ;; 
		s  ) species=$OPTARG ;;
        f  ) file=$OPTARG ;;
		i  ) index_file=$OPTARG ;;
        r  ) resource=$OPTARG ;;
		a  ) assembly=$OPTARG ;;
	    p  ) align_type='sampe' ;;
        #m  ) mask='' ;;
		h  ) echo -e $usage; return 0;;
		\? ) echo -e $usage; exit 1;;
	esac 
done

#Paramter checks
CheckVariablesOrUsage "$usage" file gender species index_home
ValidateVariableOrUsage "$usage" gender VALID_GENDERS
CheckFile $file

if [[ $file != *fastq ]]; then
	echo "Your file argument needs to be a fastq file:\t$file"
	exit
fi

lc_species=$(echo $species | tr [A-Z] [a-z])
uc_species=$(echo $species | tr [a-z] [A-Z])
index_name="${lc_species}_${gender}_${assembly}${mask}"

#Check index files
for ext in amb ann bwt pac rbwt rpac rsa sa; do
	CheckFile "${index_home}/${uc_species}/${index_name}.${ext}"
done


file_name=$(echo $file | sed 's/\.fastq$//')

echo -e "Running bwa with following options:\n\tIndex name\t= $index_name\n\tInput file\t= $file\n\tAlignment type\t= $align_type\n"

#do we need to extend the job name here to 
#make it more unique?

#This is again truncating the bsub cmd and failing?!
#bjob_cmd="bsub -J bwa_${file_name} $resource -o ${file_name}.bwa.out -e ${file_name}.bwa.err 'bwa aln $index_home/${uc_species}/${index_name}  $file > ${file_name}.${align_type}.sai; bwa $align_type ${index_home}/${uc_species}/${index_name} ${file_name}.${align_type}.sai $file > ${file_name}.${align_type}.sam'"

#Use submitJob
job_name="bwa_${file_name}"
bsub_cmd="$resource -o ${file_name}.bwa.out -e ${file_name}.bwa.err"
job_cmd="'bwa aln $index_home/${uc_species}/${index_name}  $file > ${file_name}.${align_type}.sai; bwa $align_type ${index_home}/${uc_species}/${index_name} ${file_name}.${align_type}.sai $file > ${file_name}.${align_type}.sam'"

#bjob_cmd="bsub -J bwa_${file_name} $resource -o ${file_name}.bwa.out -e ${file_name}.bwa.err 'bwa aln $index_home/${uc_species}/${index_name}  $file > ${file_name}.${align_type}.sai; bwa $align_type ${index_home}/${uc_species}/${index_name} ${file_name}.${align_type}.sai $file > ${file_name}.${align_type}.sam'"

#Add this so submitJob does not fail
#Move to efg.config?
export QUEUE_MANAGER=LSF
submitJob "$job_name" "$bsub_cmd" "$job_cmd" 


#exit

#$bjob_cmd




