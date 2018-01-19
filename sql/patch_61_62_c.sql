-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

/** 
@header patch_61_62_c.sql - feature_type.sequence_ontology
@desc   Add so_accession and so_name and fields
*/


#Or generic associated_ontology table?


ALTER table feature_type ADD `so_accession` varchar(64) DEFAULT NULL;
ALTER table feature_type ADD `so_name` varchar(255) DEFAULT NULL;

ALTER table feature_type ADD KEY `so_accession_idx` (`so_accession`);

analyze table feature_type;
optimize table feature_type;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_61_62_c.sql|feature_type.sequence_ontology');
