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
@header patch_84_85_j.sql - Normalise regulatory feature table.
@desc   Normalise regulatory feature table.
*/

DROP TABLE IF EXISTS regulatory_feature_nr;
CREATE TABLE regulatory_feature_nr (
  regulatory_feature_id int(10) unsigned NOT NULL auto_increment,
  feature_type_id int(10) unsigned default NULL,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_strand tinyint(1) NOT NULL,
  seq_region_start int(10) unsigned NOT NULL,
  seq_region_end int(10) unsigned NOT NULL,
  stable_id  varchar(16) DEFAULT NULL,
  projected boolean default FALSE,
  bound_start_length mediumint(3) unsigned NOT NULL,
  bound_end_length mediumint(3) unsigned NOT NULL,
  epigenome_count smallint(6),
  PRIMARY KEY  (regulatory_feature_id),
  UNIQUE KEY fset_seq_region_idx (seq_region_id,seq_region_start, feature_type_id),
  KEY feature_type_idx (feature_type_id),
  KEY stable_id_idx (stable_id)
) ENGINE=MyISAM;

insert into regulatory_feature_nr (
  seq_region_id,
  seq_region_start,
  seq_region_end,
  seq_region_strand,
  feature_type_id,
  stable_id,
  projected,
  bound_start_length,
  bound_end_length,
  epigenome_count 
) select 
  seq_region_id,
  seq_region_start,
  seq_region_end,
  seq_region_strand,
  feature_type_id,
  stable_id,
  projected,
  bound_start_length,
  bound_end_length,
  epigenome_count 
from regulatory_feature 
group by 
  seq_region_id,
  seq_region_start,
  seq_region_end,
  seq_region_strand,
  feature_type_id,
  stable_id,
  projected,
  bound_start_length,
  bound_end_length,
  epigenome_count 
;

DROP TABLE IF EXISTS regulatory_feature_feature_set;
create table regulatory_feature_feature_set (
  regulatory_feature_feature_set_id int(10) unsigned NOT NULL auto_increment,
  regulatory_feature_id int(10) unsigned default NULL,
  stable_id  varchar(16) DEFAULT NULL,
  feature_set_id int(10) unsigned default NULL,
  activity tinyint(3),
  PRIMARY KEY  (regulatory_feature_feature_set_id),
  UNIQUE KEY uniqueness_constraint_idx (feature_set_id,regulatory_feature_id,activity),
  KEY feature_set_idx (feature_set_id),
  KEY regulatory_feature_idx (regulatory_feature_id)
) ENGINE=MyISAM;

insert into regulatory_feature_feature_set (
  stable_id,
  feature_set_id,
  activity
) select
  stable_id,
  feature_set_id,
  activity
from regulatory_feature;

alter table regulatory_feature_feature_set add index foo (stable_id);

update regulatory_feature_feature_set, regulatory_feature_nr 
set regulatory_feature_feature_set.regulatory_feature_id=regulatory_feature_nr.regulatory_feature_id
where regulatory_feature_feature_set.stable_id=regulatory_feature_nr.stable_id;

alter table regulatory_feature_feature_set drop index foo;
alter table regulatory_feature_feature_set drop column stable_id;

drop table regulatory_feature;

rename table regulatory_feature_nr to regulatory_feature;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_j.sql|Normalise regulatory feature table.');
