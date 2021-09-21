-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
@header patch_96_97_d.sql - Fix foreign key data type inconsistencies
@desc 	Fix foreign key data type inconsistencies
*/

ALTER TABLE `chance` MODIFY COLUMN `signal_alignment_id` int(10) UNSIGNED DEFAULT NULL;
ALTER TABLE `chance` MODIFY COLUMN `control_alignment_id` int(10) UNSIGNED DEFAULT NULL;
ALTER TABLE `chance` MODIFY COLUMN `analysis_id` SMALLINT(10) unsigned DEFAULT NULL;
ALTER TABLE `segmentation` MODIFY COLUMN `regulatory_build_id` int(22) UNSIGNED DEFAULT NULL;
ALTER TABLE `motif_feature_peak` MODIFY COLUMN `motif_feature_id` int(11) unsigned NOT NULL;
ALTER TABLE `motif_feature_peak` MODIFY COLUMN `peak_id` int(11) UNSIGNED NOT NULL;
ALTER TABLE `motif_feature_regulatory_feature` MODIFY COLUMN `motif_feature_id` int(11) UNSIGNED NOT NULL;
ALTER TABLE `motif_feature_regulatory_feature` MODIFY COLUMN `regulatory_feature_id` int(11) UNSIGNED NOT NULL;
ALTER TABLE `motif_feature_regulatory_feature` MODIFY COLUMN `epigenome_id` int(11) UNSIGNED;
ALTER TABLE `binding_matrix_frequencies` MODIFY COLUMN `binding_matrix_id` int(11) UNSIGNED NOT NULL;
ALTER TABLE `binding_matrix_transcription_factor_complex` MODIFY COLUMN `binding_matrix_id` int(11) UNSIGNED NOT NULL;
ALTER TABLE `probe_set` MODIFY COLUMN `array_chip_id` int(10) UNSIGNED DEFAULT NULL;
ALTER TABLE `underlying_structure` MODIFY COLUMN `motif_feature_id` int(11) UNSIGNED NOT NULL;
ALTER TABLE `underlying_structure` MODIFY COLUMN `regulatory_feature_id` int(11) UNSIGNED NOT NULL;


-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_96_97_d.sql|Fix foreign key data type inconsistencies');
