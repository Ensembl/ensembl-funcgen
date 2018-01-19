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
@header patch_87_88_c.sql - sample_regulatory_feature_id field for regulatory build
@desc   Adds a sample_regulatory_feature_id column to the regulatory build table and puts an id in for every regulatory build.
*/

alter table regulatory_build add column sample_regulatory_feature_id int(10) unsigned;

-- Find a sample regulatory feature with many activities

drop table if exists `temp_sample_id_max_activities`;
drop table if exists `temp_sample_id`;

create table temp_sample_id_max_activities as 
select 
  regulatory_build_id, max(num_activities) as max_num_activities 
from (
  select 
    regulatory_build_id, regulatory_feature_id, count(distinct activity) num_activities 
  from 
    regulatory_build 
    join regulatory_feature using (regulatory_build_id) 
    join regulatory_activity using (regulatory_feature_id) 
  group by 
    regulatory_build_id, regulatory_feature_id 
  order by 
    num_activities desc
) a group by regulatory_build_id;

create table temp_sample_id as 
select 
  temp_sample_id_max_activities.regulatory_build_id, min(regulatory_feature_id) as sample_regulatory_feature_id 
from 
  temp_sample_id_max_activities join (
    select 
      regulatory_build_id, regulatory_feature_id, count(distinct activity) num_activities 
    from 
      regulatory_build 
      join regulatory_feature using (regulatory_build_id) 
      join regulatory_activity using (regulatory_feature_id) 
    group by 
      regulatory_build_id, regulatory_feature_id 
    order by 
      num_activities desc
) a on (
  a.regulatory_build_id=temp_sample_id_max_activities.regulatory_build_id 
  and temp_sample_id_max_activities.max_num_activities=a.num_activities
)
group by 
  temp_sample_id_max_activities.regulatory_build_id
;

update 
  regulatory_build, temp_sample_id 
set 
  regulatory_build.sample_regulatory_feature_id=temp_sample_id.sample_regulatory_feature_id 
where 
  regulatory_build.regulatory_build_id=temp_sample_id.regulatory_build_id
;

drop table `temp_sample_id_max_activities`;
drop table `temp_sample_id`;

-- Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_87_88_c.sql|sample_regulatory_feature_id field for regulatory build');
