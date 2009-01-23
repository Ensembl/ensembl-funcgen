#!/usr/local/bin/bash



type=$1
name=$2
file=$3


usage="Must provide a file_type, array_name and a file_path i.e. FASTA|TXT|CSV ARRAYNAME path/to/file"

#We could head the files first by a number of line and test the expected last line
#Then sed the header out.

if [ ! $file ]; then
	echo $usage
	exit 1;
fi

if [ $type = FASTA ]; then
	sed "s/^>/>${name}:/" $file > ${file}.tmp

elif [ $type = TXT ]; then

	#awk -F"[[:space:]]" "{print \">${name}:\" \$15 \"\n\" \$19}" $file > ${file}.tmp

	awk -F"[[:space:]]" "{if (\$1==\"Species\") probes=1;
else if (\$1==\"[Controls]\") exit;
else if (probes == 1) print \">${name}:\" \$15 \"\n\" \$19}" $file > ${file}.tmp


	#Remove spurious entries generated from header
	#sed "/^>${name}:$/d" ${file}.tmp > ${file}.tmp1
	#sed "/^$/d" ${file}.tmp1 > ${file}.tmp
	#sed "/^>${name}:Array_Address_Id$/d"  ${file}.tmp >  ${file}.tmp1 
	#sed "/^Chromosome$/d"  ${file}.tmp1 >  ${file}.tmp
	#rm -f ${file}.tmp1		

elif [ $type = CSV ]; then
			#ProbeID is 3rd field, Probe_Sequence is 10th
	awk -F"," "{print \">${name}:\" \$3 \"\n\" \$10}" $file > ${file}.tmp
	#Now remove thr first two spurious line generated from the header
	sed "/^>${name}:ProbeId$/d" ${file}.tmp > ${file}.tmp1
	sed "/^Probe_Sequence$/d" ${file}.tmp1 > ${file}.tmp
	
	rm -f ${file}.tmp1

else
	echo $usage
	exit 1;
fi

echo "head ${file}.tmp"
head ${file}.tmp

fasta_file="${file%.*}.fasta"

echo "Does this look okay? Now you need to:"
echo " mv ${file}.tmp $fasta_file"

