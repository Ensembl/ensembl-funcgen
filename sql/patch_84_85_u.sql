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
@header patch_84_85_u.sql - Remove regulatory build entries from feature_set table, relink everything else.
@desc Remove regulatory build entries from feature_set table, relink everything else.
*/

alter table regulatory_feature_feature_set add column epigenome_id int(10) unsigned default NULL;

update regulatory_feature_feature_set, feature_set
set regulatory_feature_feature_set.epigenome_id = feature_set.epigenome_id
where regulatory_feature_feature_set.feature_set_id = feature_set.feature_set_id;

alter table regulatory_feature_feature_set drop key uniqueness_constraint_idx;
alter table regulatory_feature_feature_set add unique key uniqueness_constraint_idx (epigenome_id, regulatory_feature_id);

alter table regulatory_feature_feature_set drop column feature_set_id;

rename table regulatory_feature_feature_set to regulatory_activity;

alter table regulatory_activity change regulatory_feature_feature_set_id regulatory_activity_id int(10) unsigned not null auto_increment;
alter table regulatory_evidence change regulatory_feature_feature_set_id regulatory_activity_id int(10) unsigned not null;

drop table if exists regulatory_build;
create table regulatory_build (
  regulatory_build_id int(4) unsigned not null auto_increment,
  name varchar(45)           not null,
  version                    varchar(50) DEFAULT NULL,
  initial_release_date       varchar(50) DEFAULT NULL,
  last_annotation_update     varchar(50) DEFAULT NULL,
  feature_type_id int(4)     unsigned not null,
  analysis_id int(4)         unsigned not null,
  is_current                 boolean not null default 0,
  primary key  (regulatory_build_id)
) ENGINE=MyISAM;

insert into regulatory_build (
  feature_type_id,
  analysis_id,
  version,
  is_current,
  name
) select 
  feature_type_id, 
  analysis_id, 
  "1",
  true,
  "Regulatory features" 
from feature_set where type = "regulatory" group by "feature_type_id, analysis_id";

alter table regulatory_feature add column regulatory_build_id int(10) unsigned default null;
alter table regulatory_feature add unique key uniqueness_constraint_idx (
  feature_type_id,
  seq_region_id,
  seq_region_strand,
  seq_region_start,
  seq_region_end,
  stable_id,
  bound_start_length,
  bound_end_length,
  regulatory_build_id
);

delete from feature_set where type = "regulatory";

update regulatory_feature, regulatory_build set regulatory_feature.regulatory_build_id = regulatory_build.regulatory_build_id;

insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_u.sql|Remove regulatory build entries from feature_set table, relink everything else.');
