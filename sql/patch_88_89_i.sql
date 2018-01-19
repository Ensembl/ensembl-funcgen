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
@header patch_88_89_i.sql - New columns for array table
@desc   Adds columns to array table to store data that was previously held in a module of the probemapping pipeline
*/

alter table array add column is_probeset_array       tinyint(1) NOT NULL DEFAULT 0;
alter table array add column is_linked_array         tinyint(1) NOT NULL DEFAULT 0;
alter table array add column has_sense_interrogation tinyint(1) NOT NULL DEFAULT 0;

update array set is_probeset_array       = 1 where class in ("AFFY_UTR", "AFFY_ST");
update array set is_linked_array         = 1 where class in ("AFFY_ST", "AGILENT", "PHALANX", "CODELINK", "LEIDEN", "STEMPLE_LAB_SANGER", "NIMBLEGEN_MODENCODE", "SLRI", "UCSF", "WUSTL");
update array set has_sense_interrogation = 1 where class in ("AFFY_ST");

--  Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_88_89_i.sql|New columns for array table');
