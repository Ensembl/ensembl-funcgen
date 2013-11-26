-- Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

# patch_54_55_d.sql
#
# title: edb.db_name variation fix
#
# description:
# Correct edb db_name for variation db

update external_db set db_name=replace(db_name, '_core_', '_variation_') where db_name like "%_core_Variation";



INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_d.sql|edb.db_name variation fix');

