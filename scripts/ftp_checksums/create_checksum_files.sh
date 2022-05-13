#!/bin/bash

#Function: 
#	Creates a "CHECKSUM" file recursively for each directory and subdirectories
#	containing the checksums and file names for the directory.

#Parameters: 
#	[-o]: Optional. Overwrites previously existing CHECKSUM files.
#	-p path: Mandatory. The path where the script will be creating the CHECKSUM files recursively.

function usage {
	echo "Parameters:
		[-o]: Optional. Overwrites previously existing CHECKSUM files.
		-p path: Mandatory. The path where the script will be creating the CHECKSUM files recursively."
}

function create_checksums {
  for file_name in $(ls -p $1 | grep -v /); do
    if [[ -L "$file_name" ]]; then
      continue
    fi
		if [ "$file_name" != "CHECKSUM" ]; then
		  md5sum "$file_name" >> CHECKSUM
		fi
	done
}


ini_dir=`pwd`
overwrite=0
while getopts 'op:' OPTION; do
  case "$OPTION" in
    o)
      overwrite=1;;
    p)
      path="$OPTARG";;
    ?)
      echo "script usage: $0 [-o overwrite] -p path" >&2
      usage
      exit 1;;
  esac
done

shift "$(( OPTIND - 1 ))"

if [ -z "$path" ]; then
        echo 'Path parameter is missing (-p)' >&2
        usage
	cd "$ini_dir" || exit 1
        exit 1
fi

cd "$path" || exit 1
if [ -f "CHECKSUM" ]; then
  if [ $overwrite -eq 1 ]; then
    rm "CHECKSUM"
    create_checksums "$path"
  fi
else
  create_checksums "$path"
fi

for dir in $(find "$path" -type d); do
  cd "$dir" || exit 1

  if [ -f 'CHECKSUM' ]; then
    if [ $overwrite -eq 0 ]; then
      continue
    else
      rm 'CHECKSUM'
    fi
  fi

  create_checksums "$dir"

done

cd "$ini_dir" || exit 1;

