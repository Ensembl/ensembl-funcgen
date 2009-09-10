#!/bin/sh


#Test for ENV_NAME here?

PASS=$1
shift

#Test for pass

#$* is an optional list of bed file paths
#Contents of input_dir used if not supplied


perl $EFG_SRC/scripts/import/parse_and_import.pl\
    -name  bloodomics\
    -format SEQUENCING\
	-parser Bed\
	-vendor SOLEXA\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-group efg\
	-species homo_sapiens\
	-experimental_set bloodomics\
	-ucsc_coords\
	-registry_host ensembldb.ensembl.org\
	-registry_user anonymous\
	-assembly 37\
	-host localhost\
	-port 3306\
	-dbname das_homo_sapiens_funcgen_56_37a\
	-pass $PASS\
	-cell_type CD4\
	-feature_type H3K4ac\
	-feature_analysis WindowInterval\
	-tee\
	-recover\
	-result_files $*       
