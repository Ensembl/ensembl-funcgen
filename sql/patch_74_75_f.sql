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
@header patch_74_75_f.sql - experiment.feature_cell_type_id
@desc   Added cell_type_id and feature_type_id fields to experiment table
*/


ALTER TABLE experiment ADD `feature_type_id` INT(10) unsigned         NOT NULL;
ALTER TABLE experiment ADD `cell_type_id`    INT(10) unsigned DEFAULT NULL;

ALTER TABLE experiment ADD KEY feature_type_idx(feature_type_id);
ALTER TABLE experiment ADD KEY cell_type_idx(cell_type_id);

-- Now update the id fields based on non-control input_subsets

UPDATE 
  experiment e, 
  input_subset iss 
SET 
  e.cell_type_id=iss.cell_type_id 
WHERE 
  iss.experiment_id = e.experiment_id AND 
  iss.is_control    = 0;

UPDATE 
  experiment e, 
  input_subset iss 
SET 
  e.cell_type_id    = iss.cell_type_id 
WHERE 
  iss.experiment_id = e.experiment_id AND 
  iss.is_control    = 0;

UPDATE 
  experiment e, 
  input_subset iss 
SET 
  e.feature_type_id = iss.feature_type_id 
WHERE 
  iss.experiment_id = e.experiment_id AND 
  iss.is_control    = 0;


-- sneak in a feature_set.cell_type key
ALTER TABLE feature_set ADD KEY cell_type_idx(cell_type_id);

analyze  table experiment;
optimize table experiment;
analyze  table feature_set;
optimize table feature_set;
 


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_f.sql|experiment.feature_cell_type_id');


