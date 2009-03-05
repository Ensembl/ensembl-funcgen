#!/usr/local/bin/bash


if [  -z "$EFG_SRC" ] || [ ! -d $EFG_SRC ]; then
   echo ":: You have not yet initialised the eFG environment"
   return	
fi


. $EFG_SRC/scripts/environments/funcs.sh


file_type=$1
array_name=$2
file_path=$3


usage="Usage: ARG[1] file_type ARG[2] array_name ARG[3] file_path 
i.e. pre_process_codelink.sh FASTA|TXT|CSV ARRAYNAME /path/to/file"

#We could head the files first by a number of line and test the expected last line
#Then sed the header out.

#Turn this into CheckVariablesUsage
CheckVariablesOrUsage "$usage" file_type array_name file_path



new_file_path=$(echo $file_path | sed -r 's/[[:space:]]/_/g')
new_file_path=$(echo $new_file_path | sed 's/[()]//g')
file_path=$(echo $file_path | sed 's/\([[:space:]()]\)/\\\1/g')
eval mv $file_path $new_file_path
file_path=$new_file_path
CheckFile $file_path


if [ $file_type = FASTA ]; then
	sed "s/^>/>${array_name}:/" $file_path > ${file_path}.tmp

elif [ $file_type = TXT ]; then

	#awk -F"[[:space:]]" "{print \">${name}:\" \$15 \"\n\" \$19}" $file > ${file}.tmp

	#awk "BEGIN { FS = \"\t\" }; {if (\$1==\"Species\") probes=1;
#else if (\$1==\"[Controls]\") exit;
#else if (probes == 1) print \">${name}:\" \$14 \"\n\" \$18}" $file > ${file}.tmp

awk "{if (\$1==\"CUSTOMER_PROBE_NAME\") probes=1;
else if (probes == 1) print \">${array_name}:\" \$1 \"\n\" \$2}" $file_path > ${file_path}.tmp

else
	echo $usage
	exit 1;
fi

echo "head ${file_path}.tmp"
head ${file_path}.tmp

fasta_file="${file_path%.*}.fasta"

echo "Does this look okay? Now you need to:"
echo "mv ${file_path}.tmp $fasta_file"


