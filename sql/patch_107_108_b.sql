-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
@header patch_107_108_b.sql - New epigenome_track table
@desc   New table should improve future funcgen size reduction
*/

CREATE TABLE `epigenome_track` (
  `epigenome_track_id` INT(10) unsigned NOT NULL AUTO_INCREMENT,
  `epigenome_id` INT(10) unsigned NOT NULL,
  `feature_type_id` INT(10) unsigned NOT NULL,
  `data_file_id` INT(11) unsigned NOT NULL,
  `track_type` VARCHAR(50),
  INDEX et_index ( epigenome_id, feature_type_id ),
  PRIMARY KEY (epigenome_track_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_107_108_b.sql|New epigenome_track table');