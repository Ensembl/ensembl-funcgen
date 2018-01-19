#!/bin/sh

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
ARGS=$(echo $@ | sed "s/$PASS//")



$EFG_SRC/scripts/import/import_type.pl\
	#Mandatory
	-type        "type"\                #e.g. FeatureType CellType Analysis 
    -name        "feature_name"\        #e.g. H3K36me3
    -dbname      "your_db_name"\        #e.g. your_mus_musculus_funcgen_51_37d
	-species     "latin_name"\          #e.g mus_musculus
	-description "wordy description of feature"\
    -pass $PASS\
	-class       "class of FeatureType"\ #e.g. HISTONE
	$ARGS
	#Optional
