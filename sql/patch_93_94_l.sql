-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
@header patch_93_94_l.sql - Modify binding_matrix_table
@desc Modifications related to the new motif pipeline
*/

ALTER TABLE `binding_matrix` DROP COLUMN `frequencies`;
ALTER TABLE `binding_matrix` DROP COLUMN `description`;
ALTER TABLE `binding_matrix` DROP COLUMN `feature_type_id`;
ALTER TABLE `binding_matrix` DROP COLUMN `analysis_id`;
ALTER TABLE `binding_matrix` ADD COLUMN `source` VARCHAR(20) NOT NULL;
ALTER TABLE `binding_matrix` ADD COLUMN `stable_id` VARCHAR(128) NOT NULL;
ALTER TABLE `binding_matrix` MODIFY COLUMN `name` varchar(200) NOT NULL;
ALTER TABLE `binding_matrix` DROP KEY `name_analysis_idx`;
ALTER TABLE `binding_matrix` ADD UNIQUE KEY `name_idx` (`name`);

-- patch identifier
INSERT INTO `meta` (`species_id`, `meta_key`, `meta_value`) VALUES (NULL, 'patch', 'patch_93_94_l.sql|Modify binding_matrix_table');
