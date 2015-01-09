-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
@header patch_63_64_e.sql - object_xref.ensembl_object_type 
@desc   Add Experiment to the type of ensembl_object_type field
*/

ALTER table object_xref MODIFY ensembl_object_type ENUM('Experiment', 'RegulatoryFeature', 'ExternalFeature', 'AnnotatedFeature', 'FeatureType', 'ProbeSet', 'Probe', 'ProbeFeature') not NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_e.sql|object_xref.ensembl_object_type');


