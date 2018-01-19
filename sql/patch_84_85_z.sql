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
@header patch_84_85_z.sql - Move segmentation entries from result_set table into the new segmentation_file table.
@desc Move segmentation entries from result_set table into the new segmentation_file table.
*/

drop table if exists segmentation_file;

create table segmentation_file (
  segmentation_id     int(10) unsigned NOT NULL auto_increment,
  regulatory_build_id int(10),
  name                varchar(100) default NULL,
  analysis_id         smallint(5) unsigned NOT NULL,
  epigenome_id        int(10) unsigned default NULL,
  PRIMARY KEY  (segmentation_id),
  UNIQUE KEY name_idx (name),
  KEY epigenome_idx (epigenome_id),
  KEY analysis_idx (analysis_id)
);

insert into segmentation_file (
  name,
  analysis_id,
  epigenome_id
)
select 
  name, 
  analysis_id, 
  epigenome_id 
from 
  result_set 
where 
  feature_class="segmentation";

update segmentation_file, regulatory_build set segmentation_file.regulatory_build_id=regulatory_build.regulatory_build_id;

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_z.sql|Move segmentation entries from result_set table into the new segmentation_file table.');
