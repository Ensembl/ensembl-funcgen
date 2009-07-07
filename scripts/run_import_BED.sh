#!/bin/sh


#Test for ENV_NAME here?

PASS=$1
shift

#Test for pass

#$* is an optional list of bed file paths
#Contents of input_dir used if not supplied


$EFG_SRC/scripts/import/parse_and_import.pl\
    -name  LMI_CD4_H3K4ac\
    -format SEQUENCING\
	-parser Bed\
	-vendor SOLEXA\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-group efg\
	-species homo_sapiens\
	-experimental_set LMI_CD4_H3K4ac\
	-ucsc_coords\
	-registry_host ens-staging\
	-registry_user ensro\
	-assembly 37\
	-host ens-genomics1\
	-port 3306\
	-dbname nj_importer_homo_sapiens_funcgen_55_37\
	-pass $PASS\
	-cell_type CD4\
	-feature_type H3K4ac\
	-feature_analysis WindowInterval\
	-tee\
	-recover\
	-result_files $*       
