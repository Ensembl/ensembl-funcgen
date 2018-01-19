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
@header patch_61_62_c.sql - interdb_stable_id
@desc   Add 'internal' motif_feature/external_feature stable_ids.  
        Required for inter-DB linking without using internal db ids
        i.e variation consequences
*/


ALTER table motif_feature ADD `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL;

set @n = 0;
update motif_feature mf
  join (select motif_feature_id, @n := @n + 1 new_stable_id
          from motif_feature
          order by motif_feature_id) v
    on mf.motif_feature_id = v.motif_feature_id
  set mf.interdb_stable_id = v.new_stable_id;

ALTER table motif_feature ADD UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`);


analyze table motif_feature;
optimize table motif_feature;


ALTER table external_feature ADD `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL;

set @n = 0;
update external_feature mf
  join (select external_feature_id, @n := @n + 1 new_stable_id
          from external_feature
          order by external_feature_id) v
    on mf.external_feature_id = v.external_feature_id
  set mf.interdb_stable_id = v.new_stable_id;

ALTER table external_feature ADD UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`);


analyze table external_feature;
optimize table external_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_61_62_b.sql|interdb_stable_id');
