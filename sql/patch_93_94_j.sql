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
@header patch_93_94_j.sql - Create motif_feature_peak table
@desc Stores associations between motif_features and peaks
*/

DROP TABLE IF EXISTS `motif_feature_peak`;
CREATE TABLE `motif_feature_peak` (
  `motif_feature_peak_id` int(11) NOT NULL AUTO_INCREMENT,
  `motif_feature_id` int(11) NOT NULL,
  `peak_id` int(11) NOT NULL,
  PRIMARY KEY (`motif_feature_peak_id`),
  UNIQUE KEY `motif_feature_idx` (`motif_feature_id`)
)ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO `meta` (`species_id`, `meta_key`, `meta_value`) VALUES (NULL, 'patch', 'patch_93_94_j.sql|Create motif_feature_peak table');
