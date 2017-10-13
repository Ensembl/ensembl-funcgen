-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
@header patch_motif_g.sql - Modify motif_feature table
@desc Remove display_label, rename stable_id
*/

ALTER TABLE `motif_feature` DROP COLUMN `display_label`;
ALTER TABLE `motif_feature` DROP COLUMN `interdb_stable_id`;
ALTER TABLE `motif_feature` ADD COLUMN `stable_id` VARCHAR(18) DEFAULT NULL;
ALTER TABLE `motif_feature` ADD UNIQUE KEY `unique_idx` (`binding_matrix_id`, `seq_region_id`, `seq_region_start`, `seq_region_strand`);

