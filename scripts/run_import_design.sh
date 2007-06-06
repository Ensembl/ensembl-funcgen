#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

PASS=$1


$EFG_SRC/scripts/import_design.pl\
	-location Hinxton\
	-user ensadmin\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-format TILED\
	-vendor DESIGN\
	-port 3306\
	-host ens-genomics1\
	-notes_file /lustre/work1/ensembl/nj1/efg/input/TEST_DESIGNS/DesignNotes.txt \
	-array_file /lustre/work1/ensembl/nj1/efg/input/TEST_DESIGNS/Mus_musculus.NCBIM36.40.dna.chromosome.17.+.l80.us0.t77-80.g6.p30.c148.prb \
	-dbname sanger_test_43_36e\
	-group efg\
	-data_version 42_36d\
	-tee\
	-pass $PASS\
	-recover
