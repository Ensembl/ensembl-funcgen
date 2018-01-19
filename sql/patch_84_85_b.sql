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
@header patch_84_85_b.sql - rename cell_type table to epigenome
@desc   The cell_type table may contain entries that are not specific to a single cell type only,
		i.e. brain tissue samples. Therefore the term epigenome is a more suitable one.
*/

RENAME TABLE cell_type to epigenome;
ALTER TABLE epigenome change cell_type_id epigenome_id INT(10) UNSIGNED NOT NULL auto_increment;

RENAME TABLE cell_type_lineage to epigenome_lineage;
ALTER TABLE epigenome_lineage change `cell_type_id` `epigenome_id` INT(10) UNSIGNED NOT NULL;

ALTER TABLE regulatory_feature change `cell_type_count` `epigenome_count` SMALLINT(6);

ALTER TABLE feature_set change `cell_type_id` `epigenome_id` INT(10) UNSIGNED;
ALTER TABLE feature_set drop index cell_type_idx;
ALTER TABLE feature_set add key epigenome_idx (epigenome_id);

ALTER TABLE result_set change `cell_type_id` `epigenome_id` INT(10) UNSIGNED;
ALTER TABLE result_set drop index cell_type_idx;
ALTER TABLE result_set add key epigenome_idx (epigenome_id);

ALTER TABLE input_set change `cell_type_id` `epigenome_id` INT(10) UNSIGNED;
ALTER TABLE input_set drop index cell_type_idx;
ALTER TABLE input_set add key epigenome_idx (epigenome_id);

ALTER TABLE input_subset change `cell_type_id` `epigenome_id` INT(10) UNSIGNED;

ALTER TABLE experiment change `cell_type_id` `epigenome_id` INT(10) UNSIGNED;
ALTER TABLE experiment drop index cell_type_idx;
ALTER TABLE experiment add key epigenome_idx (epigenome_id);

ALTER TABLE experimental_chip change `cell_type_id` `epigenome_id` INT(10) UNSIGNED;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_b.sql|rename cell_type table');
