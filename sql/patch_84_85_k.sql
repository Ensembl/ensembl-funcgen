-- Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
@header patch_84_85_k.sql - Normalise regulatory feature table
@desc   Link up the new non redundant regulatory features. The new regulatory_feature_ids are set. The connection is made using the stable ids.
*/

create table regulatory_feature_id_map as
select distinct 
  regulatory_feature.regulatory_feature_id as regulatory_feature_id_old, 
  regulatory_feature_nr.regulatory_feature_id  as regulatory_feature_id_new 
from 
  regulatory_feature join regulatory_feature_nr using (stable_id);

create index foo on regulatory_feature_id_map (regulatory_feature_id_old);

CREATE TABLE `regulatory_attribute_new` (
  `regulatory_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated', 'motif') default NULL,
  PRIMARY KEY  (`regulatory_feature_id`, `attribute_feature_table`, `attribute_feature_id`),
  KEY attribute_feature_idx (`attribute_feature_id`, `attribute_feature_table`)
) ENGINE=MyISAM;

insert into regulatory_attribute_new 
select 
  regulatory_feature_id_new,
  attribute_feature_id,
  attribute_feature_table
from regulatory_attribute 
join regulatory_feature_id_map on (regulatory_attribute.regulatory_feature_id=regulatory_feature_id_map.regulatory_feature_id_old)
group by 
  regulatory_feature_id_new,
  attribute_feature_id,
  attribute_feature_table;
  
drop table regulatory_feature_id_map;

insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_k.sql|Normalise regulatory feature table: Link up the new non redundant regulatory features.');
