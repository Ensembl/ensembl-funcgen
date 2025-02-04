-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
@header patch_93_94_e.sql - Create transcription_factor table
@desc Stores transcription factors and links them with feature_type
*/

DROP TABLE IF EXISTS `transcription_factor`;
CREATE TABLE `transcription_factor` (
	`transcription_factor_id` int(11) NOT NULL AUTO_INCREMENT,
	`name` varchar(120) NOT NULL,
	`feature_type_id` int(10) unsigned,
	`gene_stable_id` varchar(128),
	PRIMARY KEY (`transcription_factor_id`),
	UNIQUE KEY `name_idx` (`name`),
	KEY `feature_type_id_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO `meta` (`species_id`, `meta_key`, `meta_value`) VALUES (NULL, 'patch', 'patch_93_94_e.sql|Create transcription_factor table');
