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
@header patch_74_75_b.sql - result_set.name_unique
@desc   Patch to make name field of result_set table unique, by appending analysis
        logic_name.
*/

UPDATE 
  result_set rs, 
  analysis a 
SET
  rs.name = concat(rs.name, '_', a.logic_name)
WHERE 
  rs.analysis_id=a.analysis_id;

ALTER TABLE result_set DROP INDEX unique_idx;
ALTER TABLE result_set ADD UNIQUE KEY name_idx (name);

-- Now add back in other keys.

ALTER TABLE result_set ADD KEY cell_type_idx (cell_type_id);
ALTER TABLE result_set ADD KEY feature_type_idx (feature_type_id);
ALTER TABLE result_set ADD KEY analysis_idx (analysis_id);
ALTER TABLE result_set ADD KEY feature_class_idx (feature_class);


analyze table result_set;
optimize table result_set;

  
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_b.sql|result_set.name_unique');


