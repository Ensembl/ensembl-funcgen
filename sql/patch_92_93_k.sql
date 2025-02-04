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
@header patch_91_92_k.sql - Add table to store frip scores
@desc   Add table to store frip scores
*/

drop table if exists frip;

create table frip (
  frip_id int(10) unsigned not null auto_increment,
  peak_calling_id int(18) unsigned not null,
  frip double default null,
  total_reads int(14) default null,
  primary key (frip_id)
);

insert into frip (
  peak_calling_id,
  frip,
  total_reads
)
select 
  peak_calling_id, 
  prop_reads_in_peaks as frip, 
  total_reads 
from 
  peak_calling_qc_prop_reads_in_peaks
;

alter table frip add constraint peak_calling_id_unique unique (peak_calling_id);

drop table if exists peak_calling_qc_prop_reads_in_peaks;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_k.sql|Add table to store frip scores');
