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
@header patch_76_77_d.sql - Fix errornous feature_type_id in mirna_target_feature
@desc   A few records have an incorrect feature_type_id assinged
*/


UPDATE mirna_target_feature mtf JOIN feature_type ft ON mtf.display_label = ft.name SET mtf.feature_type_id = ft.feature_type_id;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_76_77_d.sql|Fix errornous feature_type_id in mirna_target_feature');



