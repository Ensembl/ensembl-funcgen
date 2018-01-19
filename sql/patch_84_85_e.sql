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
@header patch_84_85_e.sql - add/modify columns in input_subset table
@desc   Add/modify columns in input_subset table
*/

ALTER TABLE input_subset change replicate technical_replicate tinyint(3) unsigned DEFAULT 1 NOT NULL;
ALTER TABLE input_subset add biological_replicate tinyint(3) unsigned DEFAULT 1 NOT NULL after `name`;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_e.sql|add/modify columns in input_subset table');
