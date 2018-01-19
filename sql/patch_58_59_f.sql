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
  KEY `set_window_seq_region_idx` (`result_set_id`,`window_size`,`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1
 PARTITION BY LIST (window_size)
 (PARTITION p0 VALUES IN (0),
 PARTITION p30 VALUES IN (30),
 PARTITION p65 VALUES IN (65),
 PARTITION p130 VALUES IN (130), 
 PARTITION p260 VALUES IN (260), 
 PARTITION p450 VALUES IN (450), 
 PARTITION p648 VALUES IN (648), 
 PARTITION p950 VALUES IN (950),
 PARTITION p1296 VALUES IN (1296)
);

-- set_window_region_idx can't be UNIQUE as there may be duplicates in the 0 wsize collections
-- i.e. two or more probe features with the same start, originating from replicate probes
-- or a probe seq which is a substr of another probe.

-- Now clean status entries

DELETE s from status s, status_name sn where sn.name='RESULT_FEATURE_SET' and s.status_name_id=sn.status_name_id;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_58_59_f.sql|result_feature.partitions');


