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
@header patch_84_85_j.sql - Normalise regulatory feature table
@desc   Create a linking table between regulatory features and feature sets.
*/

--
-- Stable ids are only included here so the entries in the 
-- regulatory_feature_feature_set table can be linked up with their
-- counterparts in the regulatory_feature. The stable_ids are needed,
-- because the regulatory_feature_ids will change when the table is 
-- normalised and identical regulatory features are grouped together. 
--
-- The next patch will use the stable_ids to set the regulatory_feature_id.
-- After that the stable id column in the regulatory_feature_feature_set
-- table will be dropped.
--
create table if not exists regulatory_feature_feature_set (
  regulatory_feature_feature_set_id int(10) unsigned NOT NULL auto_increment,
  regulatory_feature_id int(10) unsigned default NULL,
  feature_set_id int(10) unsigned default NULL,
  stable_id_temp varchar(18) DEFAULT NULL,
  activity tinyint(3),
  PRIMARY KEY  (regulatory_feature_feature_set_id),
  UNIQUE KEY uniqueness_constraint_idx (feature_set_id,regulatory_feature_id),
  KEY feature_set_idx (feature_set_id),
  KEY regulatory_feature_idx (regulatory_feature_id)
) ENGINE=MyISAM;

insert into regulatory_feature_feature_set (
  stable_id_temp,
  feature_set_id,
  activity
) select
  stable_id,
  feature_set_id,
  activity
from 
  regulatory_feature
group by 
  stable_id,
  feature_set_id,
  activity
;

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_j.sql|Normalise regulatory feature table: Create a linking table between regulatory features and feature sets.');
