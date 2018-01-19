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
@header patch_84_85_m.sql - Normalise regulatory feature table
@desc   Clean up temporary columns and tables.
*/

alter table regulatory_feature_feature_set drop column stable_id_temp;

drop table regulatory_feature;
rename table regulatory_feature_nr to regulatory_feature;

drop table regulatory_attribute;
rename table regulatory_attribute_new to regulatory_attribute;

drop table regulatory_feature_id_map;

insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_m.sql|Normalise regulatory feature table: Clean up temporary columns and tables.');
