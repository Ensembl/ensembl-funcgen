#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg


#$* is list of ResultSet names to create ProbeWindowFeatures for

PASS=$1
shift


time $EFG_PERL $EFG_SRC/scripts/miscellaneous/create_transcript_result_features.pl\
	-species homo_sapiens\
	-host ens-genomics1\
    -user ensadmin\
    -exon_set\
    -no_load\
	-dbname expset_homo_sapiens_funcgen_49_36k\
	-assembly 36\
	-tee\
	-rollback\
	-pass $PASS\
	-quick\
   -slice_name chromosome:NCBI36:MT\
		$*
#	-quick\
#	-slice_name chromosome:NCBI36:MT\
