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
@header patch_84_85_i.sql - Store file types.
@desc   Store file types along with the files.
*/

alter table dbfile_registry add column file_type ENUM('BAM','BAMCOV','BIGBED','BIGWIG','VCF','CRAM');
ALTER TABLE dbfile_registry DROP PRIMARY KEY, ADD PRIMARY KEY(`table_id`, `table_name`, `file_type`);

-- We currently only have bigwig in there.
update dbfile_registry set file_type='BIGWIG';

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_i.sql|Store file types along with the files.');
