-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
@header patch_??_??_?.sql - Modify schema to accomodate new motif pipeline
@desc   
*/

CREATE TABLE `binding_matrix_frequencies` (
  `binding_matrix_frequencies_id` int(11) NOT NULL AUTO_INCREMENT,
  `binding_matrix_id` int(11) NOT NULL,
  `position` int(11) NOT NULL,
  `nucleotide` varchar(2) NOT NULL,
  `frequency` int(10) unsigned NOT NULL,
  PRIMARY KEY (`binding_matrix_frequencies_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

ALTER TABLE `binding_matrix` DROP COLUMN `frequencies`;

ALTER TABLE `binding_matrix` ADD COLUMN `source` VARCHAR(20) NOT NULL;