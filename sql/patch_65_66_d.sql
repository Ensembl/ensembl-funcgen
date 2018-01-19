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
@header patch_65_66_d.sql - Reorder an index in unmapped_object.
@desc   The unique_unmapped_obj_idx index in the unmapped_object table is
 		ineffective when querying the table.  Its first part, for example,
 		partly overlaps with the id_idx index.
*/

ALTER TABLE unmapped_object
  DROP INDEX unique_unmapped_obj_idx,
  ADD UNIQUE INDEX unique_unmapped_obj_idx
    (ensembl_id, ensembl_object_type, identifier, unmapped_reason_id,
    parent, external_db_id);

ANALYZE table unmapped_object;
OPTIMIZE table unmapped_object;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_d.sql|unmapped_object.reorder_unmapped_obj_index');
