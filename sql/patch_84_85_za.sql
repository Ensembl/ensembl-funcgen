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
@header patch_84_85_za.sql - Move entries provided by external sources from the result_set table into the new external_feature_file table.
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

-- The current names are not being used. What is displayed on the web is a 
-- concatenation of analysis_description.display_label and
-- epigenome.display_label and experimental_group.name
--
create table temp_better_names_for_external_feature_file as 
select result_set_id, concat(analysis_description.display_label, " on ", epigenome.display_label, " cells from ", experimental_group.name) as better_name_for_segmentation
from 
  result_set join experiment using (experiment_id) 
  join analysis using (analysis_id)
  join analysis_description using (analysis_id)
  join experimental_group using (experimental_group_id) 
  join epigenome on (epigenome.epigenome_id = result_set.epigenome_id)
where logic_name in (
  "WGBS_FDR_1e-4",
  "RRBS_FDR_1e-4"
);

update external_feature_file, temp_better_names_for_external_feature_file
set 
  name = temp_better_names_for_external_feature_file.better_name_for_segmentation
where
  external_feature_file.result_set_id = temp_better_names_for_external_feature_file.result_set_id;

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
drop table temp_better_names_for_external_feature_file;

-- The experiments only existed to provide a link to experimental_group to
-- build the name to display on the web. 
--
delete from experiment where experiment_id in (
  select experiment_id from result_set where feature_class="dna_methylation" 
);

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_za.sql|Move entries provided by external sources from the result_set table into the new external_feature_file table.');
