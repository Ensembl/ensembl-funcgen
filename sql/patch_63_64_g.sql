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
@header patch_63_64_g.sql - index_tidy_up
@desc   Tidy indexes on unmapped_object to reflect the core schema
*/




ALTER table unmapped_object drop KEY id_idx;
ALTER table unmapped_object drop KEY  anal_idx;
ALTER table unmapped_object drop KEY anal_exdb_idx;
ALTER table unmapped_object drop KEY object_type_idx;

ALTER table unmapped_object add UNIQUE KEY `unique_unmapped_obj_idx` (`identifier`,`ensembl_id`,`parent`,`unmapped_reason_id`,`ensembl_object_type`,`external_db_id`);
ALTER table unmapped_object add KEY `anal_exdb_idx` (`analysis_id`,`external_db_id`);
ALTER table unmapped_object add  KEY `id_idx` (`identifier`(50));

/** 
reduce identifier index to 50 characters as most ids are under 30, index will be faster
add index for queries using external_db_id in the where clause
*/

ALTER table unmapped_object add  KEY `ext_db_identifier_idx` (`external_db_id`,`identifier`);



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_g.sql|index_tidy_up');


