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
@header patch_89_90_e.sql - Table Probe: Rename array_chip_id to array_id
@desco Table Probe: Rename array_chip_id to array_id in probe and probe_set table
*/

ALTER TABLE probe change array_chip_id array_id INT(10);
ALTER TABLE probe_set change array_chip_id array_id INT(10);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_89_90_e.sql|Table Probe: Rename array_chip_id to array_id');
