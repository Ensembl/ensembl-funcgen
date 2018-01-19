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
@header patch_75_76_d.sql - feature_set.type
@desc Add mirna feature to type

*/

ALTER TABLE feature_set MODIFY type enum('annotated','regulatory','external','segmentation', 'mirna_target') DEFAULT NULL;
ALTER TABLE object_xref MODIFY ensembl_object_type enum('AnnotatedFeature','Experiment','ExternalFeature','FeatureType','MirnaTargetFeature','Probe','ProbeFeature','ProbeSet','RegulatoryFeature') NOT NULL;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_d.sql|feature_set.type mirna; object_xref.ensembl_object_type add MirnaTargetFeature');
