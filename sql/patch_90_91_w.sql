-- Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
@header patch_90_91_w.sql - 
@desc   
*/

DROP TABLE IF EXISTS `read_file_experimental_configuration`;
CREATE TABLE `read_file_experimental_configuration` (

  `read_file_experimental_configuration_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `read_file_id` int(10) unsigned,

  `experiment_id` int(10) unsigned NOT NULL,
  `biological_replicate` tinyint(3) unsigned NOT NULL DEFAULT '1',
  `technical_replicate`  tinyint(3) unsigned NOT NULL DEFAULT '1',
  
  PRIMARY KEY (`read_file_experimental_configuration_id`),
  UNIQUE KEY `name_exp_idx` (`experiment_id`, `biological_replicate`, `technical_replicate`),
  KEY `experiment_idx` (`experiment_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into read_file_experimental_configuration (
  read_file_id,
  experiment_id,
  biological_replicate,
  technical_replicate
) select 
      input_subset_id,
      experiment_id,
      biological_replicate,
      technical_replicate
from input_subset;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_w.sql|');
