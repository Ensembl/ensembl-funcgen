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
@header patch_65_66_e.sql - Add regbuild_string table
@desc   Move all meta regbuild strings to regbuild_string table 
        to avoid meta index limits on length.
*/


DROP TABLE IF EXISTS regbuild_string;

#either need to add species_id or data_set_id in here?

CREATE TABLE `regbuild_string` (
  `regbuild_string_id` int(10) NOT NULL auto_increment,
  `name` varchar(60) NOT NULL,
  `species_id` int(10) unsigned default '1',	
  `string` text NOT NULL,
  PRIMARY KEY  (`regbuild_string_id`),
  UNIQUE KEY `name_species_idx` (`species_id`, `name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- very small table with varying string sizes, so no need to set max_rows/avg_row_length for text

INSERT into regbuild_string select NULL, meta_key, species_id, meta_value from meta where meta_key like "regbuild.%_ids";

-- DELETE from meta where meta_key like "regbuild.%_ids";
-- This would make this patch unrecoverable. Do this manually

SELECT 'DELETE from meta where meta_key like "regbuild.%_ids"' as 'For safety do this part of the funcgen patch_65_66_e manually:';


# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_e.sql|add_regbuild_string_table');
