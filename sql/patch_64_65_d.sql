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
@header patch_64_65_d.sql - segmentation_feature
@desc   Add segmentation_feature table
*/


DROP TABLE IF EXISTS segmentation_feature;

CREATE TABLE `segmentation_feature` (
  `segmentation_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `feature_type_id`     int(10) unsigned default NULL,
  `feature_set_id`      int(10) unsigned default NULL,
  `score` double DEFAULT NULL,
  `display_label` varchar(60) default NULL,		
  PRIMARY KEY  (`segmentation_feature_id`),
  UNIQUE KEY `fset_seq_region_idx` (`feature_set_id`, `seq_region_id`,`seq_region_start`),
  KEY `feature_type_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;


-- score and display label are redundant here
-- display label is only ever set for external_features!

-- Is this in the base SetFeature class?

-- Add to feature_set type enum
ALTER TABLE feature_set MODIFY `type` enum('annotated','regulatory','external','segmentation') DEFAULT NULL;

-- Add to ox.ensembl_object_type enum?



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_d.sql|segmentation_feature');


