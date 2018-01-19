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

# patch_59_60_c.sql
#
# title: add motif_feature and binding_matrix
#
# description:
# Add motif_feature table, and associated link table

-- add new motif_feature table

DROP TABLE IF EXISTS `motif_feature`;
CREATE TABLE `motif_feature` (
  `motif_feature_id` int(10) unsigned NOT NULL auto_increment,
  `binding_matrix_id` INT(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) default NULL,
  `score` double default NULL,
  PRIMARY KEY  (`motif_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `binding_matrix_idx` (`binding_matrix_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- Add link table

DROP TABLE IF EXISTS `associated_motif_feature`;
CREATE TABLE `associated_motif_feature` (
   `annotated_feature_id` int(10) unsigned NOT NULL,
   `motif_feature_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`annotated_feature_id`, `motif_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



-- cell_type style query will have to be done via query extension to feature_set using cell_type_id.

DROP TABLE IF EXISTS `binding_matrix`;
CREATE  TABLE `binding_matrix` (
 `binding_matrix_id` INT(10) unsigned NOT NULL auto_increment,
 `name` VARCHAR(45) NOT NULL,
 `type` VARCHAR(45) NOT NULL,
 `feature_type_id` int(10) unsigned NOT NULL,
 `frequencies` TEXT NOT NULL,
 `description` VARCHAR(255) NULL,
 PRIMARY KEY (`binding_matrix_id`) ,
 KEY `name_type_idx` (`name`, `type`),
 KEY `feature_type_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Need to add AVG_ROW_LENGTH here for free text/blob field.
-- frequencies to be recast as blob
-- change type to analysis_id to avoid enum limitations
-- Then have wrapper methods for Jaspar/Inferred?

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_59_60_c.sql|motif_feature_binding_matrix');


