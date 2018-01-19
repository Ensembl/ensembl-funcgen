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
@header patch_84_85_l.sql - Normalise regulatory feature table
@desc   Link up the regulatory attributes with the linking table.
*/

drop table if exists regulatory_feature_id_map;

create table regulatory_feature_id_map as
select distinct 
  regulatory_feature.regulatory_feature_id as regulatory_feature_id_old, 
  regulatory_feature_nr.regulatory_feature_id  as regulatory_feature_id_new,
  regulatory_feature.feature_set_id as feature_set_id,
  regulatory_feature_feature_set.regulatory_feature_feature_set_id
from 
  regulatory_feature join regulatory_feature_nr using (stable_id)
  join regulatory_feature_feature_set on (
    regulatory_feature_feature_set.regulatory_feature_id=regulatory_feature_nr.regulatory_feature_id
    and regulatory_feature_feature_set.feature_set_id=regulatory_feature.feature_set_id
  );

create index foo on regulatory_feature_id_map (regulatory_feature_id_old);

drop table if exists regulatory_attribute_new;
  
CREATE TABLE regulatory_attribute_new (
  `regulatory_feature_feature_set_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated', 'motif') default NULL,
  PRIMARY KEY  (`regulatory_feature_feature_set_id`, `attribute_feature_table`, `attribute_feature_id`),
  KEY attribute_feature_idx (`attribute_feature_id`, `attribute_feature_table`)
) ENGINE=MyISAM;

insert into regulatory_attribute_new 
select
  regulatory_feature_feature_set_id,
  attribute_feature_id,
  attribute_feature_table
from 
  regulatory_attribute join regulatory_feature_id_map on (
    regulatory_attribute.regulatory_feature_id = regulatory_feature_id_map.regulatory_feature_id_old
  )
group by 
  regulatory_feature_feature_set_id,
  attribute_feature_id,
  attribute_feature_table
;

insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_l.sql|Normalise regulatory feature table: Link up the regulatory attributes with the linking table.');
