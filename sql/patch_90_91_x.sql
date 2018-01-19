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
@header patch_90_91_x.sql - Rename result_set to alignment in various tables and columns
@desc   Rename result_set to alignment in various tables and columns
*/

rename table result_set_input to alignment_read_file;

alter table alignment_read_file change column table_id read_file_id int(10) unsigned NOT NULL;
alter table alignment_read_file change column result_set_id alignment_id int(10) unsigned NOT NULL;
alter table alignment_read_file change column result_set_input_id alignment_read_file_id int(10) unsigned NOT NULL AUTO_INCREMENT;
alter table alignment_read_file drop column table_name;

rename table result_set to alignment;

update dbfile_registry set table_name="alignment" where table_name="result_set";

alter table alignment change column result_set_id alignment_id int(10) unsigned NOT NULL AUTO_INCREMENT;
alter table alignment drop column feature_class;
alter table alignment drop column epigenome_id;
alter table alignment drop column feature_type_id;
alter table alignment drop column experiment_id;

rename table result_set_qc_chance to alignment_qc_chance;

alter table alignment_qc_chance change column result_set_qc_chance_id alignment_qc_chance int(10) unsigned NOT NULL;
alter table alignment_qc_chance change column signal_result_set_id    signal_alignment_id int(10) unsigned NOT NULL;

rename table result_set_qc_flagstats to alignment_qc_flagstats;

alter table alignment_qc_flagstats change column result_set_qc_id alignment_qc_flagstats_id int(10) unsigned NOT NULL;
alter table alignment_qc_flagstats change column result_set_id    alignment_id              int(10) unsigned NOT NULL;

rename table result_set_qc_phantom_peak to alignment_qc_phantom_peak;

alter table alignment_qc_phantom_peak change column result_set_qc_phantom_peak_id    alignment_qc_phantom_peak_id int(10) unsigned NOT NULL;
alter table alignment_qc_phantom_peak change column result_set_id alignment_id  int(10) unsigned NOT NULL;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_x.sql|Rename result_set to alignment in various tables and columns');
