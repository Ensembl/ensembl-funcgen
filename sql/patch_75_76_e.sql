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
@header patch_75_76_e.sql - reg build adjustment
@desc add has_evidence, cell_type_count to regulatory_feature, change UNIQUE

*/
ALTER TABLE regulatory_feature ADD has_evidence TINYINT(1);
ALTER TABLE regulatory_feature ADD cell_type_count smallint(6);
ALTER TABLE regulatory_feature DROP INDEX fset_seq_region_idx;
ALTER TABLE regulatory_feature ADD UNIQUE fset_seq_region_idx (`feature_set_id`, `seq_region_id`, `seq_region_start`, `feature_type_id`);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_e.sql|add has_evidence, cell_type_count to regulatory_feature, adjust UNIQUE constraint');
