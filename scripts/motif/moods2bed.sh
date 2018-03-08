#!/usr/bin/env bash

moods_raw_dir=$1
moods_bed_dir=$moods_raw_dir/../bed

i=0
while read line
do
    mood_files[ $i ]="$line"
    (( i++ ))
done < <(ls $moods_raw_dir)

moods_raw_file=${mood_files[$LSB_JOBINDEX - 1]}
# echo ${moods_raw_file}

moods_bed_file=${moods_raw_file/.out/.bed}
# echo $moods_bed_file

awk 'BEGIN{FS=","; OFS="\t"}; { split($1, contig_info, " "); len=length($6); print contig_info[1],$3+1,$3+len,$2,"1000",$4}' $moods_raw_dir/$moods_raw_file > $moods_bed_dir/$moods_bed_file
