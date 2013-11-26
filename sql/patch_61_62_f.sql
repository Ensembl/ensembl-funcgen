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

/** 
@header patch_61_62_f.sql - regulatory_feature.fset_seq_region_idx
@desc   Add a unique key to the regulatory_feature table and remove 
        other unecessary keys
*/


ALTER TABLE regulatory_feature ADD unique KEY `fset_seq_region_idx` (`feature_set_id`, `seq_region_id`,`seq_region_start`);

-- Drop other keys as there are no Slice methods which do not use feature_set_ids
ALTER TABLE regulatory_feature DROP KEY `seq_region_idx`;
ALTER TABLE regulatory_feature DROP KEY `feature_set_idx`;

optimize table regulatory_feature;
analyze table regulatory_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_61_62_f.sql|regulatory_feature.fset_seq_region_idx');


