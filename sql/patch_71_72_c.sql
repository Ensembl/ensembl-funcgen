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
@header patch_71_72_c.sql - supporting_set PK
@desc   Add 'type' to PK of supporting_set PK
*/

ALTER TABLE supporting_set DROP PRIMARY KEY, ADD PRIMARY KEY(`data_set_id`,`supporting_set_id`,`type`);

OPTIMIZE TABLE supporting_set;
ANALYZE TABLE supporting_set;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
 VALUES (NULL, 'patch', 'patch_71_72_c.sql|added_type_to_supporting_set_PK');
