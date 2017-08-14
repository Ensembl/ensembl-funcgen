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
@header patch_90_91_i.sql - Translate sequence region ids of annotated features
@desc   Translate sequence region ids of annotated features
*/

DROP TABLE IF EXISTS `temp_annotated_feature_with_translated_ids`;

CREATE TABLE `temp_annotated_feature_with_translated_ids` (
  `annotated_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  `score` double DEFAULT NULL,
  `feature_set_id` int(10) unsigned NOT NULL,
  `summit` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`annotated_feature_id`),
  UNIQUE KEY `seq_region_feature_set_idx` (`seq_region_id`,`seq_region_start`,`feature_set_id`),
  KEY `feature_set_idx` (`feature_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=39;

insert into temp_annotated_feature_with_translated_ids
select
	annotated_feature_id,
	core_seq_region_id as seq_region_id,
	seq_region_start,
	seq_region_end,
	seq_region_strand,
	display_label,
	score,
	feature_set_id,
	summit
from annotated_feature join seq_region using (seq_region_id);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_i.sql|Translate sequence region ids of annotated features');
