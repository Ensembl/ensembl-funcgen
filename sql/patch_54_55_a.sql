# patch_54_55_a.sql
#
# title: Partition result_feature table
#
# description:
# Partition result-feature table based on a window_size, seq_region_id key
# Initial arbitrary partition number of 250, should see most seq_region_id/window_size 
# keys hosted in their own partition (244 for current human).

--DROP TABLE IF EXISTS `part_result_feature`;
CREATE TABLE `part_result_feature` (
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
   PARTITION BY KEY (`window_size`,`seq_region_id`)
   PARTITIONS 250;



INSERT into part_result_feature (SELECT * from result_feature);

DROP TABLE IF EXISTS `result_feature`;
CREATE TABLE `part_result_feature` (
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
   PARTITION BY KEY (`window_size`,`seq_region_id`)
   PARTITIONS 250;


INSERT into result_feature (SELECT * from part_result_feature);

DROP TABLE `part_result_feature`;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_a.sql|partition_result_feature');


