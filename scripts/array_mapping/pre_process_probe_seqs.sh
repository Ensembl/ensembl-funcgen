#!/usr/local/bin/bash

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


if [  -z "$EFG_SRC" ] || [ ! -d $EFG_SRC ]; then
   echo ":: You have not yet initialised the eFG environment"
   return	
fi


. $EFG_SRC/scripts/environments/funcs.sh

#opt this as we don't need vendor for FASTA format
vendor=$1
file_type=$2
array_name=$3
file_path=$4

export VALID_VENDORS='ILLUMINA CODELINK PHALANX AGILENT AFFY'
#AGILENT is 

usage="Usage: ARG[1] vendor ARG[2] file_type ARG[3] array_name ARG[4] file_path 
i.e. pre_process_seqs.sh ILLUMINA|CODELINK|PHALANX|AGILENT FASTA|TXT|CSV|GEO ARRAYNAME /path/to/file"

#We could head the files first by a number of line and test the expected last line
#Then sed the header out.

vendor=$(echo $vendor | tr [:lower:] [:upper:])
file_type=$(echo $file_type | tr [:lower:] [:upper:])
CheckVariablesOrUsage "$usage" file_type array_name file_path


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
	ValidateVariableOrUsage "$usage" vendor VALID_VENDORS

	#awk -F"[[:space:]]" "{print \">${name}:\" \$15 \"\n\" \$19}" $file > ${file}.tmp

	if [ $vendor = ILLUMINA ]; then
		awk "BEGIN { FS = \"\t\" }; {if (\$1==\"Species\") probes=1;
else if (\$1==\"[Controls]\") exit;
else if (probes == 1) print \">${array_name}:\" \$14 \"\n\" \$18}" $file_path > ${file_path}.tmp

	elif [ $vendor = CODELINK ]; then
awk "{if (\$1==\"CUSTOMER_PROBE_NAME\") probes=1;
else if (probes == 1) print \">${array_name}:\" \$1 \"\n\" \$2}" $file_path > ${file_path}.tmp

	elif [ $vendor = PHALANX ]; then
		awk "{if (\$1==\"Phalanx_PrbInx\") probes=1;
else if (probes == 1) print \">${array_name}:\" \$1 \"\n\" \$2}" $file_path > ${file_path}.tmp

  elif [ $vendor = AFFY ]; then
     awk "{if (\$1==\"Probe\") probes=1;
else if (probes == 1) print \">probe:${array_name}:\" \$1 \":\" \$2 \":\" \$3 \";\" \"\n\" \$5}" $file_path > ${file_path}.tmp
      
	else
		echo "No support for $vendor $file_type file"
		exit 1
	fi

elif [ $file_type = CSV ]; then

	if [ $vendor = ILLUMINA ]; then
			#ProbeID is 3rd field, Probe_Sequence is 10th
		#echo '	awk -F"," "{print \">${name}:\" \$3 \"\n\" \$10}" '

		#awk -F"," "{print \">${name}:\" \$3 \"\n\" \$10}" $file > ${file}.tmp



		awk "BEGIN { FS = \",\" }; {if (\$1==\"Search_key\") probes=1;
#else if (probes == 1) print \">${array_name}:\" \$3 \"\n\" \$10}" $file_path > ${file_path}.tmp




		#Methylation arrays headers
		#IllmnID or Name & SourceSeq
		#Nice standardised formats?
		
		#Just hard code for now

		#echo "hardcoded for infinium 27"

		#Infinium i.e. HumanMethylation450
		#awk "BEGIN { FS = \",\" }; { if (\$1==\"IlmnID\") probes=1;
#else if (probes == 1) print \">${array_name}:\" \$1 \"\n\" \$14 }" $file_path > ${file_path}.tmp


#awk "BEGIN { FS = \",\" }; { if (\$1==\"Name\") probes=1;
#else if (probes == 1) print \">${array_name}:\" \$1 \"\n\" \$7 }" $file_path > ${file_path}.tmp



	#Now remove thr first two spurious line generated from the header
		#sed "/^>${name}:ProbeId$/d" ${file}.tmp > ${file}.tmp1
		#sed "/^Probe_Sequence$/d" ${file}.tmp1 > ${file}.tmp
	
		#rm -f ${file}.tmp1
	else
		echo "No support for $vendor $file_type file, maybe you want FASTA file type?"
		exit 1
	fi
elif [ $file_type = GEO ]; then

    awk "BEGIN { FS = \"\t\" }; {if (\$1==\"ID\") probes=1; else if ((probes == 1) && (\$6==\"FALSE\")) print \">probe:${array_name}:\" \$4 \";\" \"\n\" \$20}" $file_path > ${file_path}.tmp


else
	echo $usage
	exit 1;
fi


#Delete blank lines
sed '/^[:space:]*$/d' ${file_path}.tmp > ${file_path}.temp
rm -f ${file_path}.tmp

echo "head ${file_path}.temp"
head ${file_path}.temp

fasta_file="${file_path%.*}.fasta"

echo "
Does this look okay? Now you need to:"
echo "mv ${file_path}.temp $fasta_file"


