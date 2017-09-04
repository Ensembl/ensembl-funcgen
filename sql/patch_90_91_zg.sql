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
@header patch_90_91_zg.sql - 
@desc   
*/

rename table dbfile_registry to data_file;

ALTER TABLE data_file ADD data_file_id INT NOT NULL AUTO_INCREMENT unique first;

ALTER TABLE data_file DROP PRIMARY KEY;
ALTER TABLE data_file add PRIMARY KEY (data_file_id);

alter table alignment add column bam_file_id    int DEFAULT null;
alter table alignment add column bigwig_file_id int DEFAULT null;

update 
  alignment, data_file 
set 
  alignment.bam_file_id=data_file_id 
where 
  file_type="bam" 
  and table_id=alignment_id 
  and table_name="alignment"
;

update 
  alignment, data_file 
set 
  alignment.bigwig_file_id=data_file_id 
where 
  file_type="bigwig" 
  and table_id=alignment_id 
  and table_name="alignment"
;

-- Check like this:
-- select * from alignment join data_file on (data_file.table_name="alignment" and data_file.table_id=alignment.alignment_id) where bam_file_id != data_file_id and file_type="BAM"

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_zg.sql|');
