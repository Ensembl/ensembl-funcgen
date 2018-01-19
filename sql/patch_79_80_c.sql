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
@header patch_79_80_c.sql - stable_id_change_to_varchar
@desc   Changed to regulatory_feature.stable_id field to a varchar from an int
*/

ALTER table regulatory_feature MODIFY `stable_id` varchar(128) DEFAULT NULL;

# data patch for those species which have stable_ids
# update regulatory_feature set stable_id=lpad(stable_id, 11,'0');
# update regulatory_feature set stable_id=concat('ENSR', stable_id); # Or ENSMUSR for mouse

SELECT "WARNING: Need to apply species specific stable_id prefix patch manually";
optimize table regulatory_feature;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_79_80_c.sql|stable_id_changed_to_varchar');

