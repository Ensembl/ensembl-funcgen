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

# patch_54_55_a.sql
#
# title: Partition result_feature table
#
# description:
# Partition result-feature table based on a window_size, seq_region_id key
# Initial arbitrary partition number of 250, should see most seq_region_id/window_size 
# keys hosted in their own partition (244 for current human).
# To do this on a species by species basis just use those seq_region_ids which are 
# available in the probe_feature table for this release
# This will mean that any subsequent seq_regions/windows are added after 
# partitioning, will share a partition.

--DROP TABLE IF EXISTS `part_result_feature`;
#CREATE TABLE `part_result_feature` (
#  `result_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
#  `result_set_id` int(10) unsigned NOT NULL,
#  `seq_region_id` int(10) unsigned NOT NULL,
#  `seq_region_start` int(10) NOT NULL,
#  `seq_region_end` int(10) NOT NULL,
#  `seq_region_strand` tinyint(4) NOT NULL,
#  `window_size` smallint(5) unsigned NOT NULL,
#  `score` double DEFAULT NULL,
#	KEY (result_feature_id),
#   KEY `set_window_seq_region_idx` (`result_set_id`,`window_size`,`seq_region_id`,`seq_region_start`)
#   ) ENGINE=MyISAM DEFAULT CHARSET=latin1 
#   PARTITION BY KEY (`window_size`,`seq_region_id`)
#   PARTITIONS 250;

#INSERT into part_result_feature (SELECT * from result_feature);

rename table result_feature to tmp_result_feature;



DROP TABLE IF EXISTS `result_feature`;
CREATE TABLE `result_feature` (
  `result_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `result_set_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) NOT NULL,
  `seq_region_end` int(10) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `window_size` smallint(5) unsigned NOT NULL,
  `score` double DEFAULT NULL,
	KEY (result_feature_id),
   KEY `set_window_seq_region_idx` (`result_set_id`,`window_size`,`seq_region_id`,`seq_region_start`)
   ) ENGINE=MyISAM DEFAULT CHARSET=latin1 
   PARTITION BY KEY (`window_size`)
   PARTITIONS 7;


#For human we have 93 seq_regions in result_feature
#7 windows gives 644 partitions, way too many
#Just stay with num windows for DBs with no data
#For DBs with data do 100, or product of seq_region_ids * num windows
#InnoDB cannot be used due to copying problems


INSERT into result_feature (SELECT * from tmp_result_feature);

DROP TABLE `tmp_result_feature`;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_a.sql|partition_result_feature');


