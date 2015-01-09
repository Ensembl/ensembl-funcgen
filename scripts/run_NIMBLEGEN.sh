#!/bin/sh

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

