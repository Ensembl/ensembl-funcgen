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
@header patch_85_86_d.sql - Add new columns to input_subset table to accommodate paired-end data
@desc   Add new columns to input_subset table to accommodate paired-end data
*/

ALTER TABLE input_subset ADD COLUMN `read_length` int(10) DEFAULT NULL;
ALTER TABLE input_subset ADD COLUMN `is_paired_end` tinyint(1) DEFAULT NULL;
ALTER TABLE input_subset ADD COLUMN `paired_with` int(10) DEFAULT NULL;
ALTER TABLE input_subset ADD COLUMN `file_size` bigint(20) DEFAULT NULL;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_85_86_d.sql|Add new columns to input_subset table to accommodate paired-end data');