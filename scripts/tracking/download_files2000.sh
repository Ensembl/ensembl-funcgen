TOTLINES=$(wc -l < input_file_xtra_filtered_human_only_known_tagets_rel_94_8.tsv)
TOTLINES=$((TOTLINES-1))
echo "Total Lines: "$TOTLINES

count=0
while IFS=$'\t' read -r -a myArray
do
 LINK=$(echo "${myArray[24]}")
 FILE=$(echo "${myArray[11]}")

 if [[ $count -gt 2000 && $count -lt 3001 ]]
 then
 	if [ ! -e $FILE ]
 	then
	    echo "Downloading file: " $LINK
	    wget $LINK -P /nfs/production/panda/ensembl/funcgen/source_files/fastq/homo_sapiens/encode/v94/
	    echo "Downloaded " $count " of " $TOTLINES
	    #echo $LINK
 	fi
 
     #break
 fi
 count=$((count+1))
 
done < input_file_xtra_filtered_human_only_known_tagets_rel_94_8.tsv 


