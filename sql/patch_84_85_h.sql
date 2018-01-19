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
@header patch_84_85_h.sql - Store file types.
@desc   Store file types along with the files.
*/

alter table dbfile_registry add column file_type ENUM('BAM','BAMCOV','BIGBED','BIGWIG','VCF','CRAM', 'DIR');
alter table dbfile_registry drop primary key, add primary key(table_id, table_name, file_type);
alter table dbfile_registry add column md5sum varchar(45) default null;

-- Some file names have trailing whitespaces, hence using RTRIM
update dbfile_registry set file_type='BIGWIG' where RIGHT(RTRIM(path),3) = ".bw";
update dbfile_registry set file_type='BIGBED' where RIGHT(RTRIM(path),3) = ".bb";


-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_h.sql|Store file types along with the files.');
