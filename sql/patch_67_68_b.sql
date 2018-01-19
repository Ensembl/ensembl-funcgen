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
@header patch_67_68_b.sql - input_subset.archive_id_display_url
@desc   Move the archive_id and data_url fields from experiment to input_subset
*/


-- Move currently misplaced fields from experiment table
ALTER TABLE input_subset ADD `archive_id` varchar(20) DEFAULT NULL;
ALTER TABLE input_subset ADD `display_url` varchar(255) DEFAULT NULL;

-- Migrate data from experiment
UPDATE input_subset iss, input_set iset, experiment e set iss.archive_id = e.archive_id where e.experiment_id=iset.experiment_id and iset.input_set_id=iss.input_set_id;
UPDATE input_subset iss, input_set iset, experiment e set iss.display_url = e.data_url where e.experiment_id=iset.experiment_id and iset.input_set_id=iss.input_set_id;

-- Probably need to do some clean here wrt to pooled reps, as these will all have been assigned the same archive_id
-- Should have others in tracking DB

ALTER TABLE input_subset ADD INDEX `archive_idx`(`archive_id`);

OPTIMIZE TABLE input_subset;
ANALYZE TABLE input_subset;

-- Need to drop these fields from experiment
ALTER TABLE experiment DROP archive_id; -- also drops index
ALTER TABLE experiment DROP data_url;

OPTIMIZE TABLE experiment;
ANALYZE TABLE experiment;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_67_68_b.sql|input_subset.archive_id_display_url');


