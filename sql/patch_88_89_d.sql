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

/**
@header patch_88_89_d.sql - create probe_transcript table
@desc   Creates probe_transcript table, moves data from object_xref and xref into it.
*/

DROP TABLE IF EXISTS `probe_transcript`;

CREATE TABLE `probe_transcript` (
  `probe_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_id`    int(10) unsigned NOT NULL,
  `stable_id`   varchar(128)      NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_transcript_id`),
  KEY `probe_transcript_id` (`probe_transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into `probe_transcript` (`probe_id`, `stable_id`, `description`) (
  select 
    ensembl_id as probe_id, 
    dbprimary_acc as stable_id, 
    linkage_annotation as description 
  from 
    object_xref 
    join xref using (xref_id) 
  where ensembl_object_type="Probe"
);

--  Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_88_89_d.sql|created probe_transcript table');
