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
@header patch_90_91_r.sql - Translate coord_system_ids in meta_coord table
@desc   Translate coord_system_ids in meta_coord table
*/

delete from coord_system where coord_system.is_current=false;

DROP TABLE IF EXISTS `temp_meta_coord_with_translated_ids`;

CREATE TABLE `temp_meta_coord_with_translated_ids` (
  `table_name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `max_length` int(11) DEFAULT NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into temp_meta_coord_with_translated_ids 
select 
  meta_coord.table_name,
  coord_system.core_coord_system_id,
  min(meta_coord.max_length)
from meta_coord join coord_system using (coord_system_id)
group by 
  meta_coord.table_name, coord_system.core_coord_system_id
order by 
  meta_coord.table_name
;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_r.sql|Translate coord_system_ids in meta_coord table');
