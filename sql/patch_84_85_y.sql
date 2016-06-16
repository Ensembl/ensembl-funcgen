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
@header patch_84_85_y.sql - Move entries provided by external sources from the result_set table into the new external_feature_file table.
@desc Move entries provided by external sources from the result_set table into the new external_feature_file table.
*/

drop table if exists external_feature_file;

create table external_feature_file (
  external_feature_file_id  int(10) unsigned NOT NULL auto_increment,
  name                      varchar(100) default NULL,
  analysis_id               smallint(5) unsigned NOT NULL,
  epigenome_id              int(10) unsigned default NULL,
  feature_type_id           int(10) unsigned default NULL,
  experiment_id             int(10) unsigned default NULL,
  result_set_id             int(10) unsigned default NULL,
  PRIMARY KEY  (external_feature_file_id),
  UNIQUE KEY name_idx (name),
  KEY epigenome_idx (epigenome_id),
  KEY analysis_idx (analysis_id)
);

insert into external_feature_file (
  name,
  analysis_id,
  epigenome_id,
  feature_type_id,
  experiment_id,
  result_set_id
)
select 
  name, 
  analysis_id, 
  epigenome_id, 
  feature_type_id,
  experiment_id,
  result_set_id
from 
  result_set 
where feature_class="dna_methylation";

-- Adjust the table_names and table_ids in the dbfile_registry table
--
update 
  dbfile_registry, external_feature_file 
set 
  dbfile_registry.table_name="external_feature_file", 
  dbfile_registry.table_id=external_feature_file.external_feature_file_id 
where 
  dbfile_registry.table_id=external_feature_file.result_set_id 
  and dbfile_registry.table_name="result_set";

delete from result_set where feature_class="dna_methylation";

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_y.sql|Move entries provided by external sources from the result_set table into the new external_feature_file table.');
