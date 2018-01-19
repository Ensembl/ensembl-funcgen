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

-- rename oligo tables and dependent oligo columns
-- drop design_type.experiment_id
 


-- oligo_probe > probe
CREATE TABLE `probe` (
  `probe_id` int(10) unsigned NOT NULL auto_increment,
   `probe_set_id` int(10) unsigned default NULL,
   `name` varchar(40) NOT NULL default '',
   `length` smallint(6) unsigned NOT NULL default '0',
   `array_chip_id` int(10) unsigned NOT NULL default '0',
   `class` varchar(20) default NULL,
    PRIMARY KEY  (`probe_id`, `name`),
    KEY `probe_set_idx` (`probe_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

INSERT into probe SELECT * from oligo_probe;
DROP table oligo_probe;

-- oligo_probe_set > probe_set
CREATE TABLE `probe_set` (
   `probe_set_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(20) NOT NULL default '',
   `size` smallint(6) unsigned NOT NULL default '0',
   `family` varchar(20) default NULL,
   PRIMARY KEY  (`probe_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

INSERT into probe_set SELECT * from oligo_probe_set;
DROP table oligo_probe_set;



-- oligo_feature > probe_feature
CREATE TABLE `probe_feature` (
   `probe_feature_id` int(10) unsigned NOT NULL auto_increment,
   `seq_region_id` int(10) unsigned NOT NULL default '0',
   `seq_region_start` int(10) NOT NULL default '0',
   `seq_region_end` int(10) NOT NULL default '0',
   `seq_region_strand` tinyint(4) NOT NULL default '0', 
   `coord_system_id` int(10) unsigned NOT NULL default '0',
   `probe_id` int(10) unsigned NOT NULL default '0',
   `analysis_id` int(10) unsigned NOT NULL default '0',	
   `mismatches` tinyint(4) NOT NULL default '0',
   `cigar_line` text,
   PRIMARY KEY  (`probe_feature_id`),
   KEY `probe_idx` (`probe_id`),
   KEY `seq_region_idx` (`seq_region_id`, `seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


INSERT into probe_feature SELECT * from oligo_feature;
DROP table oligo_feature;

-- result fields
DROP INDEX oligo_probe_idx on result;
ALTER table result CHANGE oligo_probe_id probe_id int(10) unsigned NOT NULL default '0';
CREATE INDEX probe_idx on result (`probe_id`);



-- change design_type tables

CREATE TABLE `experimental_design_type` (
   `design_type_id` int(10) unsigned NOT NULL auto_increment,
   `table_name` varchar(40) default NULL,
   `table_id` int(10) unsigned default NULL,	
   PRIMARY KEY  (`design_type_id`, `table_name`, `table_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- would need to migrate entris if there were any

alter table design_type drop column experiment_id;
