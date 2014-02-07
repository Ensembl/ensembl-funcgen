#!/bin/sh

# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#. ~/src/ensembl-efg/scripts/.efg

PASS=$1
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/parse_and_import.pl\
	-name  H3K4me1-GM06990-new\
	-format tiled\
	-location Hinxton\
	-contact njohnson@ebi.ac.uk\
	-species homo_sapiens\
	-fasta\
	-port 3306\
	-host ens-genomics1\
	-dbname sanger_test_43_36e\
	-array_name "ENCODE3.1.1"\
	-group efg\
	-vendor sanger\
	-assembly 36\
	-exp_date 2006-11-02\
	-tee\
	-verbose\
	-pass $1\
	-recover\
	$ARGS
