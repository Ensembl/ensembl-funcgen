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
@header patch_68_69_c.sql - regbuild_string.species_id_not_null
@desc   Change the regbuild_string.species_id to be not null smallint
*/

-- Note species_id is simply an int in meta, but probably needs changing

ALTER table regbuild_string MODIFY species_id smallint(5) unsigned NOT NULL DEFAULT 1;

OPTIMIZE TABLE regbuild_string;
ANALYZE TABLE regbuild_string;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_68_69_c.sql|regbuild_string.species_id_not_null');


