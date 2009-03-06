#!/usr/local/bin/bash


if [  -z "$EFG_SRC" ] || [ ! -d $EFG_SRC ]; then
   echo ":: You have not yet initialised the eFG environment"
   return	
fi


. $EFG_SRC/scripts/environments/funcs.sh


vendor=$1
file_type=$2
array_name=$3
file_path=$4

export VALID_VENDORS='ILLUMINA CODELINK PHALANX'

usage="Usage: ARG[1] vendor ARG[2] file_type ARG[3] array_name ARG[4] file_path 
i.e. pre_process_codelink.sh ILLUMINA|CODELINK|PHALANX FASTA|TXT|CSV ARRAYNAME /path/to/file"

#We could head the files first by a number of line and test the expected last line
#Then sed the header out.

vendor=$(echo $vendor | tr [:lower:] [:upper:])
file_type=$(echo $file_type | tr [:lower:] [:upper:])
CheckVariablesOrUsage "$usage" file_type array_name file_path
ValidateVariableOrUsage "$usage" vendor VALID_VENDORS


new_file_path=$(echo "$file_path" | grep -E "[[:space:]()]")

if [ "$new_file_path" ]; then
	new_file_path=$(echo $file_path | sed -r 's/[[:space:]]/_/g')
	new_file_path=$(echo $new_file_path | sed 's/[()]//g')
	file_path=$(echo $file_path | sed 's/\([[:space:]()]\)/\\\1/g')
	CheckFile "$file_path"
	eval mv $file_path $new_file_path
	file_path=$new_file_path
else
	CheckFile $file_path
fi


if [ $file_type = FASTA ]; then
	echo "Prefixing $vendor fasta header with $array_name"
	sed "s/^>/>${array_name}:/" $file_path > ${file_path}.tmp

elif [ $file_type = TXT ]; then

	#awk -F"[[:space:]]" "{print \">${name}:\" \$15 \"\n\" \$19}" $file > ${file}.tmp

	if [ $vendor = ILLUMINA ]; then
		awk "BEGIN { FS = \"\t\" }; {if (\$1==\"Species\") probes=1;
else if (\$1==\"[Controls]\") exit;
else if (probes == 1) print \">${name}:\" \$14 \"\n\" \$18}" $file > ${file}.tmp

	elif [ $vendor = CODELINK ]; then
awk "{if (\$1==\"CUSTOMER_PROBE_NAME\") probes=1;
else if (probes == 1) print \">${array_name}:\" \$1 \"\n\" \$2}" $file_path > ${file_path}.tmp

	elif [ $vendor = PHALANX ]; then
		awk "{if (\$1==\"Phalanx_PrbInx\") probes=1;
else if (probes == 1) print \">${array_name}:\" \$1 \"\n\" \$2}" $file_path > ${file_path}.tmp

	else
		echo "No support for $vendor $file_type file"
		exit 1
	fi

elif [ $type = CSV ]; then

	if [ $vendor = ILLUMINA ]; then
			#ProbeID is 3rd field, Probe_Sequence is 10th
		awk -F"," "{print \">${name}:\" \$3 \"\n\" \$10}" $file > ${file}.tmp
	#Now remove thr first two spurious line generated from the header
		sed "/^>${name}:ProbeId$/d" ${file}.tmp > ${file}.tmp1
		sed "/^Probe_Sequence$/d" ${file}.tmp1 > ${file}.tmp
	
		rm -f ${file}.tmp1
	else
		echo "No support for $vendor $file_type file"
		exit 1
	fi
else
	echo $usage
	exit 1;
fi

echo "head ${file_path}.tmp"
head ${file_path}.tmp

fasta_file="${file_path%.*}.fasta"

echo "Does this look okay? Now you need to:"
echo "mv ${file_path}.tmp $fasta_file"


