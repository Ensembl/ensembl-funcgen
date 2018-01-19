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
@header patch_74_75_c.sql - input_subset.analysis_id_experiment_idx 
@desc   An analysis_id has been added to the input_subset table.
        An experiment_id index has been added.
*/


-- This will mirror the input_set.analysis_id. Consequently, InputSubset has been 
-- changed to inherit from Set, and the feature_type validation in of Set subclass 
-- constructors has been moved to the Set constructor.
-- This work is a prerequisite to the retiring the InputSet class/table.

ALTER TABLE input_subset ADD `analysis_id` smallint(5) unsigned NOT NULL;
ALTER TABLE input_subset ADD KEY analysis_idx (analysis_id);
ALTER TABLE input_subset ADD KEY experiment_idx (experiment_id);


UPDATE 
  input_subset iss, 
  input_set_input_subset isiss, 
  input_set inp 
SET 
  iss.analysis_id  = inp.analysis_id 
WHERE 
  inp.input_set_id      = isiss.input_set_id AND
  isiss.input_subset_id = iss.input_subset_id;

analyze table input_subset;
optimize table input_subset;

  
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_c.sql|input_subset.analysis_id_experiment_idx');


