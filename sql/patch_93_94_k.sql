-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
@header patch_93_94_k.sql - Create motif_feature_regulatory_feature table
@desc Stores associations between motif features and regulatory features
*/

DROP TABLE IF EXISTS `motif_feature_regulatory_feature`;
CREATE TABLE `motif_feature_regulatory_feature` (
  `motif_feature_regulatory_feature_id` int(11) NOT NULL AUTO_INCREMENT,
  `motif_feature_id` int(11) NOT NULL,
  `regulatory_feature_id` int(11) NOT NULL,
  `epigenome_id` int(11),
  `has_matching_Peak` tinyint(3) unsigned DEFAULT '0',
  PRIMARY KEY (`motif_feature_regulatory_feature_id`),
  UNIQUE KEY `mf_rf_ep_idx` (`motif_feature_id`,`regulatory_feature_id`, `epigenome_id`)
)ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO `meta` (`species_id`, `meta_key`, `meta_value`) VALUES (NULL, 'patch', 'patch_93_94_k.sql|Create motif_feature_regulatory_feature table');
