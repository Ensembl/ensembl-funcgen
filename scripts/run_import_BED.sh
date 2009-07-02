#!/bin/sh


#Test for ENV_NAME here?

#$* is list of bed filepaths

PASS=$1
shift


$EFG_SRC/scripts/import/parse_and_import.pl\
	-name  LMI_CD4_H3K4ac\                        #experiment/input_dir name
	-format SEQUENCING\                           #experimental format
	-parser Bed\                                  #Import parser type
	-vendor SOLEXA\                               #experiment technology vendor(Will use this parser if -parser not set)
	-location Hinxton\                            #group location
	-contact njohnson@ebi.ac.uk\                  #group contact
	-group efg\                                   #group name
	-species homo_sapiens\
	-experimental_set LMI_CD4_H3K4ac\              
	-ucsc_coords\                                 #Starts at 0 rather than 1
	-registry_host ens-staging\                   #Registry params for selecting
	-registry_user ensro\                         #correct dnadb
	-assembly 37\                                 #Genome assembly version
	-host ens-genomics1\                          #DB connection params
	-port 3306\
	-dbname nj_importer_homo_sapiens_funcgen_55_37\
	-pass $PASS\	
	-cell_type CD4\                               #Feature/Cell/Analysis info
	-feature_type H3K4ac\
	-feature_analysis WindowInterval\
	-tee\                                         #Tees all log ouput to STDOUT(logfile also written)
	-recover\                                     #Required for Nimblegen 2nd stage or if experiment already exists
	$*                                            #List of file bed file paths
