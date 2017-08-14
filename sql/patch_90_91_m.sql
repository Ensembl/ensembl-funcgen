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
@header patch_90_91_m.sql - Translate sequence region ids of mi rna target features
@desc   Translate sequence region ids of mi rna target features
*/

DROP TABLE IF EXISTS `temp_mirna_target_feature_with_translated_ids`;

CREATE TABLE `temp_mirna_target_feature_with_translated_ids` (
  `mirna_target_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `feature_set_id` int(10) unsigned NOT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `accession` varchar(60) DEFAULT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  `evidence` varchar(60) DEFAULT NULL,
  `interdb_stable_id` int(10) unsigned DEFAULT NULL,
  `method` varchar(60) DEFAULT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `supporting_information` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`mirna_target_feature_id`),
  UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

insert into temp_mirna_target_feature_with_translated_ids 
select
	mirna_target_feature_id,
	feature_set_id,
	feature_type_id,
	accession,
	display_label,
	evidence,
	interdb_stable_id,
	method,
	core_seq_region_id as seq_region_id,
	seq_region_start,
	seq_region_end,
	seq_region_strand,
	supporting_information
from mirna_target_feature join seq_region using (seq_region_id);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_m.sql|Translate sequence region ids of mi rna target features');
