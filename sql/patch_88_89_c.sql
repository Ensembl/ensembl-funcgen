-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
@header patch_88_89_c.sql - create probe_feature_transcript table
@desc   Creates probe_feature_transcript table, moves data from object_xref and xref into it.
*/

drop table if exists `temp_probe_feature_transcript`;

create table temp_probe_feature_transcript as 
select 
  ensembl_id as probe_feature_id, 
  max(version) as max_version
from 
  object_xref 
  join xref using (xref_id) 
where 
  ensembl_object_type="ProbeFeature"
group by 
  ensembl_id;

drop table if exists `probe_feature_transcript`;

CREATE TABLE `probe_feature_transcript` (
  `probe_feature_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_feature_id` int(10) unsigned DEFAULT NULL,
  `stable_id`   varchar(128) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_feature_transcript_id`),
  KEY `probe_feature_transcript_id_idx` (`probe_feature_transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


insert into `probe_feature_transcript` (`probe_feature_id`, `stable_id`, `description`) (
select 
  ensembl_id as probe_feature_id, 
  dbprimary_acc as stable_id, 
  linkage_annotation as description
from 
  temp_probe_feature_transcript
  join object_xref on (probe_feature_id=ensembl_id)
  join xref on (object_xref.xref_id=xref.xref_id and temp_probe_feature_transcript.max_version=version) 
where 
  ensembl_object_type="ProbeFeature"
);

-- create unique index probe_feature_id_idx on probe_feature_transcript(probe_feature_id);
create index probe_feature_id_idx on probe_feature_transcript(probe_feature_id);

drop table if exists `temp_probe_feature_transcript`;

--  Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_88_89_c.sql|created probe_feature_transcript table');
