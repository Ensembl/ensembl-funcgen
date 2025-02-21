-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
@header patch_99_100_e.sql - Make some foreign keys not mandatory
@desc   Make some foreign keys not mandatory
*/

ALTER TABLE segmentation_cell_tables MODIFY control_alignment_id int(23) unsigned DEFAULT NULL;
UPDATE segmentation_cell_tables SET control_alignment_id = NULL WHERE control_alignment_id = 0;

ALTER TABLE object_xref MODIFY analysis_id smallint(5) unsigned DEFAULT NULL;
UPDATE object_xref SET analysis_id = NULL WHERE analysis_id = 0;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_99_100_e.sql|Make some foreign keys not mandatory');
