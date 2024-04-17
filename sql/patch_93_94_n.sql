-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
@header patch_93_94_n.sql - Modify indices in motif_feature_peak table
@desc Remove incorrect index, add two new ones
*/

ALTER TABLE `motif_feature_peak` DROP INDEX `motif_feature_idx`;
ALTER TABLE `motif_feature_peak` ADD INDEX `motif_feature_idx` (`motif_feature_id`);
ALTER TABLE `motif_feature_peak` ADD INDEX `peak_idx` (`peak_id`);

-- patch identifier
INSERT INTO `meta` (`species_id`, `meta_key`, `meta_value`) VALUES (NULL, 'patch', 'patch_93_94_n.sql|Modify indices in motif_feature_peak table');
