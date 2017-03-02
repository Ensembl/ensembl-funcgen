-- Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
@header patch_88_89_b.sql - Create probe_seq table
@desc   Creates a table for storing probe sequences.
*/

--
-- Table structure for table `probe_seq`
--

DROP TABLE IF EXISTS `probe_seq`;

CREATE TABLE `probe_seq` (
  `probe_seq_id` int(10) NOT NULL AUTO_INCREMENT,
  `probe_sha1` char(40) NOT NULL,
  `probe_dna` text NOT NULL,
  PRIMARY KEY (`probe_seq_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

ALTER TABLE probe ADD COLUMN probe_seq_id int(10) DEFAULT NULL;

--  Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_88_89_b.sql|seq_region_name_255');
