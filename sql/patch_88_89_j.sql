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
@header patch_88_89_j.sql - Added array_chip_id column to probe_set table
@desc   A column to distinguish between probe sets with the same name but on different array chips
*/

ALTER TABLE probe ADD COLUMN probe_seq_id int(10) DEFAULT NULL;
alter table probe add index `probe_seq_idx` (`probe_seq_id`); 

--  Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_88_89_j.sql|Added array_chip_id column to probe_set table');
