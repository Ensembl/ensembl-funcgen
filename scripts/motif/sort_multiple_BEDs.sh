#!/usr/bin/env bash

unsorted_bed_dir=$1
sorted_bed_dir=$unsorted_bed_dir/../sorted_bed

# echo $unsorted_bed_dir
# echo $sorted_bed_dir

# unsorted_bed_files=()
# for file in unsorted_bed_dir/*.bed
# do
#     unsorted_bed_files=("${unsorted_bed_files[@]}" "$file")
# done

i=0
while read line
do
    unsorted_bed_files[ $i ]="$line"        
    (( i++ ))
done < <(ls $unsorted_bed_dir)



unsorted_bed_file=${unsorted_bed_files[$LSB_JOBINDEX - 1]}
# echo $unsorted_bed_file
sorted_bed_file=${unsorted_bed_file/.bed/.sorted.bed}
# echo $sorted_bed_file

bedtools sort -i $unsorted_bed_dir/$unsorted_bed_file >$sorted_bed_dir/$sorted_bed_file
