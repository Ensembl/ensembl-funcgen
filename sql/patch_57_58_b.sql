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

# patch_57_58_b.sql
#
# title: result_feature_partitions
#
# description:
# Change the partition definition to list to avoid uneven modulus derived
# partition assignment,
# Also change window sizes, this may break current data as the old wsizes 
# may not match those listed in the partitions. Hence may nee to reimport
# result_feature collections.


DROP table if exists list_result_feature;

CREATE TABLE `list_result_feature` (
  `result_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `result_set_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) NOT NULL,
  `seq_region_end` int(10) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `window_size` smallint(5) unsigned NOT NULL,
  `scores` longblob NOT NULL,
  KEY `result_feature_idx` (`result_feature_id`),
  KEY `set_window_seq_region_idx` (`result_set_id`,`window_size`,`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM  DEFAULT CHARSET=latin1
PARTITION BY LIST (window_size)
 (PARTITION p0 VALUES IN (0),
 PARTITION p50 VALUES IN (50),
 PARTITION p150 VALUES IN (150), 
 PARTITION p300 VALUES IN (300), 
 PARTITION p450 VALUES IN (450), 
 PARTITION p600 VALUES IN (600), 
 PARTITION p750 VALUES IN (750),
 PARTITION p900 VALUES IN (900),	
 PARTITION p1150 VALUES IN (1150)
);


-- Will we get failures with the large blobs here?

INSERT into list_result_feature select * from result_feature;

DROP table result_feature;

RENAME table list_result_feature to result_feature;

	
-- 0 to capture natural resolution for tiling arrays (omited for seq data)
-- 50 to for high res seq displays (could theoretically be 30bp with default server/client max_allowed_packet_size?)
-- 150 intervals to capture 2-3 50bp probes/seq in one bin  
-- Max range is roughly 1Mb in 900 (best fit for 15" screen) ~ 1111 > 1150


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_57_58_b.sql|result_feature_paritions');


