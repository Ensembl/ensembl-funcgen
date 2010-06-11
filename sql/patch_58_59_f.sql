# patch_58_59_f.sql
#
# title: result_feature.partitions
#
# description:
# Alter the partitions to match the new window sizes 
# This will destroy all existing data


DROP table if exists result_feature;

CREATE TABLE `result_feature` (
  `result_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `result_set_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) NOT NULL,
  `seq_region_end` int(10) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `window_size` smallint(5) unsigned NOT NULL,
  `scores` longblob NOT NULL,
  KEY `result_feature_idx` (`result_feature_id`),
  UNIQUE KEY `set_window_seq_region_idx` (`result_set_id`,`window_size`,`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1
 PARTITION BY LIST (window_size)
 (PARTITION p0 VALUES IN (0),
 PARTITION p50 VALUES IN (30),
 PARTITION p150 VALUES IN (130), 
 PARTITION p300 VALUES IN (260), 
 PARTITION p450 VALUES IN (450), 
 PARTITION p600 VALUES IN (648), 
 PARTITION p750 VALUES IN (950),
 PARTITION p900 VALUES IN (1296)
);

-- Also Added UNIQUE key here

-- Now clean status entries

DELETE s from status s, status_name sn where sn.name='RESULT_FEATURE_SET' and s.status_name_id=sn.status_name_id;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_58_59_e.sql|result_feature.partitions');


