-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
@header patch_91_92_f.sql - Updates to alignment table
@desc   Updates to alignment table
*/

alter table alignment add column experiment_id int(15) unsigned default null;
alter table alignment add column has_duplicates boolean default null;
alter table alignment add column is_control     boolean default null;
alter table alignment add column source_alignment_id int(22) unsigned default NULL;
alter table alignment add column deduplicated_alignment_id int(28) unsigned default NULL;
alter table alignment add column to_gender enum('male','female','hermaphrodite','mixed','unknown') DEFAULT NULL;
alter table alignment add column is_complete boolean default null;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_f.sql|Updates to alignment table');
