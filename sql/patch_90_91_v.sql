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
@header patch_90_91_v.sql - Create read_file table and populate it
@desc   
*/

DROP TABLE IF EXISTS `read_file`;
CREATE TABLE `read_file` (

  `read_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(300) NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `is_paired_end` tinyint(1) DEFAULT NULL,
  `paired_with` int(10) DEFAULT NULL,
  `file_size` bigint(20) DEFAULT NULL,
  `read_length` int(10) DEFAULT NULL,
  `md5sum` varchar(45) DEFAULT NULL,
  `file` text,
  `notes` text,
  PRIMARY KEY (`read_file_id`),
  UNIQUE KEY `read_file_id_idx` (`read_file_id`)
  
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into read_file (
  read_file_id,
  name,
  analysis_id,
  is_paired_end,
  paired_with,
  file_size,
  read_length
) select 
      input_subset_id,
      name,
      analysis_id,
      is_paired_end,
      paired_with,
      file_size,
      read_length
from input_subset;

-- When patching a tracking database, use this command instead of the above:
-- 
-- insert into read_file (
--   read_file_experimental_configuration_id,
--   name,
--   analysis_id,
--   is_paired_end,
--   paired_with,
--   file_size,
--   read_length,
--   md5sum,
--   file,
--   notes
-- ) select 
--       input_subset_id,
--       name,
--       analysis_id,
--       is_paired_end,
--       paired_with,
--       file_size,
--       read_length,
--       md5sum,
--       local_url,
--       notes
-- from input_subset join input_subset_tracking using (input_subset_id);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_v.sql|Create read_file table and populate it');
