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

# patch_59_60_h.sql
#
# Title: af_amf.index_tweaks
#   
#
# Description:
#   Alter annotated_feature and motif_feature indexes to increase performance


ALTER TABLE associated_motif_feature ADD KEY `motif_feature_idx` (`motif_feature_id`);
OPTIMIZE TABLE associated_motif_feature;

ALTER TABLE annotated_feature ADD UNIQUE KEY `seq_region_feature_set_idx` (`seq_region_id`,`seq_region_start`,`feature_set_id`);
ALTER TABLE annotated_feature DROP KEY seq_region_idx; 
OPTIMIZE TABLE annotated_feature;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_59_60_h.sql|af_amf.index_tweaks');
