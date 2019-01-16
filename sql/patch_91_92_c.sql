-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2019] EMBL-European Bioinformatics Institute

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
@header patch_91_92_c.sql - Create underlying_structure table
@desc   This table associates regulatory features to motif features.
*/

CREATE TABLE `underlying_structure` (
  `underlying_structure_id` int(11) NOT NULL AUTO_INCREMENT,
  `regulatory_feature_id` int(11) NOT NULL,
  `motif_feature_id` int(11) NOT NULL,
  PRIMARY KEY (`underlying_structure_id`)
)ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO `meta` (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_91_92_c.sql|Create underlying_structure table');
