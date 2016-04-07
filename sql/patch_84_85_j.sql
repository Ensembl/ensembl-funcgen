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
@header patch_84_85_j.sql - Make activity an enum.
@desc   Make activity an enum.
*/

alter table regulatory_feature add column activity_as_enum ENUM('INACTIVE', 'REPRESSED', 'POISED', 'ACTIVE', 'NA');

update regulatory_feature set activity_as_enum = 
case activity
    when 0 then 'INACTIVE'
    when 1 then 'ACTIVE'
    when 2 then 'POISED'
    when 3 then 'REPRESSED'
    when 4 then 'NA'
    else null
end;

alter table regulatory_feature drop column activity;
alter table regulatory_feature change activity_as_enum activity ENUM('INACTIVE', 'REPRESSED', 'POISED', 'ACTIVE', 'NA');

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_j.sql|Make activity an enum.');
