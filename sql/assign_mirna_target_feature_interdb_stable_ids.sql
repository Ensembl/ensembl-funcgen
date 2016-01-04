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
@header Assign_new_interdb_stable_ids
@desc   Creates new interdb_stable_id values for NULL
        values, starting at the next logical value.
*/


-- No point in doing this as we can't guarantee they won't be reused between releases
-- if we selectively delete some of the last IDs
-- set @n =(select max(interdb_stable_id) +1 from motif_feature);

update mirna_target_feature set interdb_stable_id=NULL;
set @n = 0;

update mirna_target_feature mf
  join (select mirna_target_feature_id, @n := @n + 1 new_stable_id
          from mirna_target_feature where interdb_stable_id is NULL
          order by mirna_target_feature_id) v
    on mf.mirna_target_feature_id = v.mirna_target_feature_id
  set mf.interdb_stable_id = v.new_stable_id;

