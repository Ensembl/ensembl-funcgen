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

# patch_58_59_d.sql
#
# title: add regulatory_feature.binary_string_projected
#
# description:
# Add binary_string and projected field and also redefine table to remove 
# patch binary_string and remove unecessary AVG_ROW_LENGTH


DROP TABLE IF EXISTS `new_regulatory_feature`;
CREATE TABLE `new_regulatory_feature` (
  `regulatory_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(80) default NULL,
  `feature_type_id`	int(10) unsigned default NULL,
  `feature_set_id`	int(10) unsigned default NULL,
  `stable_id` mediumint(8) unsigned default NULL,
  `bound_seq_region_start` int(10) unsigned NOT NULL,	
  `bound_seq_region_end` int(10) unsigned NOT NULL,
  `binary_string` varchar(255) default NULL,
  `projected` boolean default FALSE,	
  PRIMARY KEY  (`regulatory_feature_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `stable_id_idx` (`stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


INSERT into new_regulatory_feature select *, NULL, FALSE from regulatory_feature;
UPDATE new_regulatory_feature set binary_string=display_label;
UPDATE new_regulatory_feature set display_label=NULL;

DROP table regulatory_feature;
RENAME table new_regulatory_feature TO regulatory_feature;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_58_59_d.sql|regulatory_feature.binary_string_projected');


