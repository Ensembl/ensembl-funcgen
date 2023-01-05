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
@header patch_93_94_m.sql - Create binding_matrix_frequencies table
@desc Stores the frequency values of a binding matrix
*/

DROP TABLE IF EXISTS `binding_matrix_frequencies`;
CREATE TABLE `binding_matrix_frequencies` (
  `binding_matrix_frequencies_id` int(11) NOT NULL AUTO_INCREMENT,
  `binding_matrix_id` int(11) NOT NULL,
  `position` int(11) unsigned NOT NULL,
  `nucleotide` enum('A','C','G','T') NOT NULL,
  `frequency` int(10) unsigned NOT NULL,
  PRIMARY KEY (`binding_matrix_frequencies_id`),
  KEY `binding_matrix_id_idx` (`binding_matrix_id`),
  UNIQUE KEY `unique_constraint_idx` (`binding_matrix_id`,`position`,`nucleotide`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO `meta` (`species_id`, `meta_key`, `meta_value`) VALUES (NULL, 'patch', 'patch_93_94_m.sql|Create binding_matrix_frequencies table');
