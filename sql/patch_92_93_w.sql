-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
@header patch_92_93_w.sql - Modify index name_exp_idx from table read_file_experimental_configuration
@desc   Modify index name_exp_idx from table read_file_experimental_configuration
*/

ALTER table read_file_experimental_configuration DROP INDEX name_exp_idx;
ALTER table read_file_experimental_configuration ADD UNIQUE INDEX name_exp_idx (experiment_id, biological_replicate, technical_replicate, paired_end_tag, multiple);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_w.sql|Modify index name_exp_idx from table read_file_experimental_configuration');
