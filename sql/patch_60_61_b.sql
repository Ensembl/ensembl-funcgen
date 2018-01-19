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

# patch_60_61_b.sql
#
# title:  add binding_matrix.analysis_id and drop binding_matrix.type
#
# description:
# Add analysis_id field to binding_matrix table and drop type field


ALTER table binding_matrix DROP KEY `name_type_idx`;

ALTER table binding_matrix ADD `analysis_id` int(10) unsigned NOT NULL;
ALTER table binding_matrix DROP `type`;

ALTER table binding_matrix ADD KEY `name_analysis_idx` (`name`, `analysis_id`);

OPTIMIZE table binding_matrix;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_60_61_b.sql|binding_matrix.analysis_id');


