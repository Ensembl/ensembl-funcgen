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

# patch_52_53_a.sql
#
# title: Unmapped object/reason tables
#
# description:
# Add and alter the unmapped tables for eFG

---Don't drop first just in case one already exists

CREATE TABLE `unmapped_object` (
  `unmapped_object_id` int(10) unsigned NOT NULL auto_increment,
  `type` enum('xref', 'probe2transcript') NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `external_db_id` smallint(5) unsigned default NULL,
  `identifier` varchar(255) NOT NULL,
  `unmapped_reason_id` smallint(5) unsigned NOT NULL,
  `query_score` double default NULL,
  `target_score` double default NULL,
  `ensembl_id` int(10) unsigned default '0',
  `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType', 'Probe', 'ProbeSet') NOT NULL,
  `parent` varchar(255) default NULL,
  PRIMARY KEY  (`unmapped_object_id`),
  KEY `id_idx` (`identifier`),
  KEY `anal_idx` (`analysis_id`),
  KEY `anal_exdb_idx` (`analysis_id`,`external_db_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


CREATE TABLE `unmapped_reason` (
  `unmapped_reason_id` smallint(5) unsigned NOT NULL auto_increment,
  `summary_description` varchar(255) default NULL,
  `full_description` varchar(255) default NULL,
  PRIMARY KEY  (`unmapped_reason_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_a.sql|unmapped_tables');


