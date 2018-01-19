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

$EFG_SRC/scripts/import_array_from_fasta.pl \
		-format CUSTOM \
		-vendor NIMBLEGEN \
		-array ARRAY_NAME \
		-dbname species_name_funcgen_74_ASSEMBLYVERSION \
    -pass $1 \
    -file Array.fasta \
		-outdir /out/dir/ \
		$ARGS
