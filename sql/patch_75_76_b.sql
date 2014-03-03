-- Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
@header patch_75_76_b.sql - result_set/experiment.display_url/archive_id
@desc   Move input_subset display_url and archive_id to the experiment table, and 
        add an experiment_id field to result_set. Drop redundant input_subset fields.
*/


ALTER TABLE experiment ADD COLUMN archive_id varchar(60) DEFAULT NULL; 
-- To handle concated archive_ids

ALTER TABLE experiment ADD COLUMN display_url varchar(255) DEFAULT NULL;
ALTER TABLE result_set ADD COLUMN experiment_id int(10) unsigned default NULL;
ALTER TABLE result_set ADD KEY    experiment_idx (experiment_id);

-- This update will likely be over-written by the input_subset patch script
-- to ensure faithful migration of all data from the old tracking schema

UPDATE experiment e, input_subset iss, result_set_input rsi, result_set rs 
  SET e.display_url=iss.display_url, e.archive_id=iss.archive_id, rs.experiment_id=e.experiment_id
  WHERE e.experiment_id=iss.experiment_id and iss.is_control=0 and iss.input_subset_id=rsi.table_id
  AND   rsi.table_name='input_subset' and rsi.result_set_id=rs.result_set_id;

analyze table result_set;
optimize table result_set;


ALTER TABLE input_subset DROP COLUMN archive_id;
-- we could maintain this, but we are currently putting the SRR ids in the name field

ALTER TABLE input_subset DROP COLUMN display_url;

  
-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_b.sql|result_set/experiment.display_url/archive_id');


