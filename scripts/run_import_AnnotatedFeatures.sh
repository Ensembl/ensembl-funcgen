#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

PASS=$1
shift

input_files=$@

#Now check is defined
#Importer currently only handles 1 file/import

if [[]]; then
fi






# Input set defines a distinct set of data used during the import. This enables import of various 
# experiments without having to integrate the data model (tiling array support is no longer maintianed). 
# This also allows references to and import tracking of data which is not necessarily imported 
# into the DB i.e. raw reads

# The cell_type, feature_type and analysis must all be predefined in the DB.
# See scripts/run_import_type.sh. Analysis entries are currently imported manually
# or via the CreateDB function which imports scripts/imports/typee/Analysis.txt


name=MyExperiment1             #Name of experiment
format=SEQUENCING              #input_set format
group=efg                      #Name of experimental_group
location=Hinxton               #Location of experimental_group
contact=njohnson@ebi.ac.uk     #contact for experimental_group
species=homo_sapiens
input_set=MyExperiment1        #Name of the input_set
input_set_feature_class=annotated  #Type of feature
reg_dbhost=reg_mysql_host      #Can alternatively defined dnadb params directly
reg_dbuser=reg_mysql_user      #params directly. See below
dbhost=mysql_host
dbname=my_homo_sapiens_funcgen_65_37
dbpot=3306
ctype=CD4
ftype=H3K4ac
analysis=MyPeakAnalysis
parser=Bed

#Optional
#vendor=' -vendor Illumina '   # Defaults to parser format
#assembly=' -assembly 37 '   # Defaults to current
ucsc_coords=' -ucsc_coords ' # accounts for half open format
tee=' -tee '                 # tees all log output to STDOUT
recover=' -recover '         # ???
#dnadb_host
#dnadb_name
#dnadb_port
#dnadb_user
#dnadb_pass


perl $EFG_SRC/scripts/import/parse_and_import.pl\
	-name             $name\
	-format           $format\
	-location         $location\
	-contact          $contact\
	-group            $group\
	-species          $species\
	-input_set        $input_set\
	-registry_host    $reg_dbhost\
	-registry_user    $reg_dbuser\
	-port             $dbport\
	-host             $dbhost\
	-dbname           $dbname\
	-cell_type        $ctype\
	-feature_type     $ftype\
	-feature_analysis $analysis\
	-parser           $parser\
	-pass             $PASS\
	$vendor\
	$assembly\
	$ucsc_coords\
	$tee\
	$recover\
	$input_files

