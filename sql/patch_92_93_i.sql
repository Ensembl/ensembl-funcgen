-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
@header patch_91_92_i.sql - Add table to store fastqc outcomes
@desc   Add table to store fastqc outcomes
*/

drop table if exists fastqc;
create table fastqc (
  fastqc_id    int(11) unsigned NOT NULL AUTO_INCREMENT,
  read_file_id int(14) unsigned NOT NULL,
  basic_statistics             enum('PASS','WARN','FAIL') default null,
  per_base_sequence_quality    enum('PASS','WARN','FAIL') default null,
  per_tile_sequence_quality    enum('PASS','WARN','FAIL') default null,
  per_sequence_quality_scores  enum('PASS','WARN','FAIL') default null,
  per_base_sequence_content    enum('PASS','WARN','FAIL') default null,
  per_sequence_gc_content      enum('PASS','WARN','FAIL') default null,
  per_base_n_content           enum('PASS','WARN','FAIL') default null,
  sequence_length_distribution enum('PASS','WARN','FAIL') default null,
  sequence_duplication_levels  enum('PASS','WARN','FAIL') default null,
  overrepresented_sequences    enum('PASS','WARN','FAIL') default null,
  adapter_content              enum('PASS','WARN','FAIL') default null,
  kmer_content                 enum('PASS','WARN','FAIL') default null,
  PRIMARY KEY (`fastqc_id`)
);

ALTER TABLE fastqc ADD CONSTRAINT read_file_id_unique UNIQUE (read_file_id);

-- insert into fastqc (
--   read_file_id,
--   basic_statistics,
--   per_base_sequence_quality,
--   per_tile_sequence_quality,
--   per_sequence_quality_scores,
--   per_base_sequence_content,
--   per_sequence_gc_content,
--   per_base_n_content,
--   sequence_length_distribution,
--   sequence_duplication_levels,
--   overrepresented_sequences,
--   adapter_content,
--   kmer_content
-- ) 
--   select 
--     distinct
--     read_file_fastqc.read_file_id,
--     basic_statistics.status,
--     per_base_sequence_quality.status,
--     per_tile_sequence_quality.status,
--     per_sequence_quality_scores.status,
--     per_base_sequence_content.status,
--     per_sequence_gc_content.status,
--     per_base_n_content.status,
--     sequence_length_distribution.status,
--     sequence_duplication_levels.status,
--     overrepresented_sequences.status,
--     adapter_content.status,
--     kmer_content.status
--   from 
--     read_file_fastqc
--     join read_file_fastqc as basic_statistics             using (read_file_id)
--     join read_file_fastqc as per_base_sequence_quality    using (read_file_id)
--     join read_file_fastqc as per_tile_sequence_quality    using (read_file_id)
--     join read_file_fastqc as per_sequence_quality_scores  using (read_file_id)
--     join read_file_fastqc as per_base_sequence_content    using (read_file_id)
--     join read_file_fastqc as per_sequence_gc_content      using (read_file_id)
--     join read_file_fastqc as per_base_n_content           using (read_file_id)
--     join read_file_fastqc as sequence_length_distribution using (read_file_id)
--     join read_file_fastqc as sequence_duplication_levels  using (read_file_id)
--     join read_file_fastqc as overrepresented_sequences    using (read_file_id)
--     join read_file_fastqc as adapter_content              using (read_file_id)
--     join read_file_fastqc as kmer_content                 using (read_file_id)
--   where
--     basic_statistics.title                  = "Basic Statistics"
--     and  per_base_sequence_quality.title    = "Per base sequence quality"
--     and  per_tile_sequence_quality.title    = "Per tile sequence quality"
--     and  per_sequence_quality_scores.title  = "Per sequence quality scores"
--     and  per_base_sequence_content.title    = "Per base sequence content"
--     and  per_sequence_gc_content.title      = "Per sequence GC content"
--     and  per_base_n_content.title           = "Per base N content"
--     and  sequence_length_distribution.title = "Sequence Length Distribution"
--     and  sequence_duplication_levels.title  = "Sequence Duplication Levels"
--     and  overrepresented_sequences.title    = "Overrepresented sequences"
--     and  adapter_content.title              = "Adapter Content"
--     and  kmer_content.title                 = "Kmer Content";

alter table fastqc add column run_failed boolean default false;
alter table fastqc add column error_message text;

drop table if exists read_file_fastqc;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_i.sql|Add table to store fastqc outcomes');
