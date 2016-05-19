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
@header patch_84_85_v.sql - Created production names for feature types.
@desc Created production names for feature types, histone names are used for generating bigwig files.
*/

ALTER TABLE feature_type add production_name VARCHAR(120) after description;

-- Create production names from the names.
--
update feature_type set production_name=name;

update feature_type set production_name = replace(production_name, ':', '_');
update feature_type set production_name = replace(production_name, '+', '');
update feature_type set production_name = replace(production_name, '-', '_');
update feature_type set production_name = replace(production_name, '.', '_');
update feature_type set production_name = replace(production_name, '/', '_');
update feature_type set production_name = replace(production_name, ' ', '_');

insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_v.sql|Created production names for feature types, histone names are used for generating bigwig files.');
