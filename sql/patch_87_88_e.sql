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
@header patch_87_88_e.sql - create probeset_transcript table
@desc   Creates probeset_transcript table, moves data from object_xref and xref into it.
*/

DROP TABLE IF EXISTS `probeset_transcript`;

CREATE TABLE `probeset_transcript` (
  `probeset_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probeset_id`    int(10) unsigned NOT NULL,
  `stable_id`   varchar(18)      NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probeset_transcript_id`),
  KEY `probeset_transcript_id_idx` (`probeset_transcript_id`)
);

insert into `probeset_transcript` (`probeset_id`, `stable_id`, `description`) (
  select 
    ensembl_id as probeset_id, 
    dbprimary_acc as stable_id, 
    linkage_annotation as description 
  from 
    object_xref 
    join xref using (xref_id) 
  where ensembl_object_type="ProbeSet"
);



