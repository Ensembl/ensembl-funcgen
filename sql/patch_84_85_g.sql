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
@header patch_84_85_g.sql - update external_db_id data type
@desc   Update external_db_id data type to align regulation with core schema
*/

alter table external_db modify external_db_id INTEGER UNSIGNED NOT NULL AUTO_INCREMENT;
alter table unmapped_object modify external_db_id INTEGER UNSIGNED;
alter table xref modify external_db_id INTEGER UNSIGNED;


-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_g.sql|update external_db_id data type');
