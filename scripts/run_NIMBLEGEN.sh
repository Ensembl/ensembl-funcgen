#!/bin/sh

PASS=$1
shift

perl $EFG_SRC/scripts/import/parse_and_import.pl\
	-name 'DVD_OR_EXPERIMENT_NAME'\
	-format tiled\
	-location Hinxton\
	-contact 'your_email@address.com'\
	-vendor NIMBLEGEN\
	-species homo_sapiens\
	-fasta\
	-port 3306\
	-host dbhost\
	-dbname 'dbname_funcgen_VERSION_BUILD'\
	-array_set\
	-array_name 'DESIGN_NAME'\
	-group efg\
	-assembly 36\
	-tee\
	-pass $PASS\
	-recover

