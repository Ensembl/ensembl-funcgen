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
@header patch_84_85_c.sql - new columns to the epigenome table
@desc   Modify existing or add new columns to the epigenome table
*/

ALTER TABLE epigenome change efo_id ontology_accession VARCHAR(20);
ALTER TABLE epigenome add ontology ENUM('EFO','CL') after ontology_accession;
ALTER TABLE epigenome add production_name VARCHAR(120) after description;

-- Create production names from the names.
--
update epigenome set production_name=name;

update epigenome set production_name = replace(production_name, ':', '_');
update epigenome set production_name = replace(production_name, '+', '');
update epigenome set production_name = replace(production_name, '-', '_');
update epigenome set production_name = replace(production_name, '.', '_');
update epigenome set production_name = replace(production_name, '/', '_');
update epigenome set production_name = replace(production_name, ' ', '_');

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_c.sql|new epigenome table columns');
