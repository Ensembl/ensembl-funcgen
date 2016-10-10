#!/bin/sh

# Using the script from here:
# https://github.com/EMBL-EBI-GCA/gca-tools/blob/master/submissions/file-manifest

usage="Usage: $(basename "$0") [-h | --help] [path]...
Generate a file manifest for submitting to GCA at EMBL-EBI
-h | --help display this help and exit
This script generates a file manifest for files that you want to submit to GCA
at EMBL-EBI. The script takes a list of directories or files to add to the
manifest. The default directory is your current working directory.
We recommend you gather the files you want to submit into a single directory,
and then write the output from this script to a file so that you can send it to
GCA at EMBL-EBI.
./file-manifest /path/to/directory > manifest.tsv
Please read your manifest before you send it to EMBL-EBI to confirm that it
contains only the files you wish to send to us.
This script is a wrapper for the following unix command. You may prefer to copy
and paste this command as an alternative to running this executable:
find -type f -printf '%p\t%s\t' -execdir sh -c 'md5sum "{}" | sed s/\ .*//' \;
"
while getopts ":h-:" opt; do
case $opt in
-) case "$OPTARG" in
help) echo "$usage"
exit
;;
*) printf "invalid option: --%s\n" "$OPTARG" >&2
echo "$usage"
exit 1
;;
esac ;;
h) echo "$usage"
exit
;;
\?) printf "invalid option: -%s\n" "$OPTARG" >&2
echo "$usage"
exit 1
;;
esac
done
shift $(($OPTIND - 1))
find $@ -xtype f -printf '%p\t%s\t' -execdir sh -c 'md5sum "{}" | sed s/\ .*//' \;