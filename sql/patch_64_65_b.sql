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
@header patch_64_65_b.sql - feature_type.analysis_id
@desc   Add analysis_id field and add to unique key. To support SegmentationFeature states
*/

ALTER table feature_type ADD `analysis_id` smallint(5) unsigned default NULL AFTER class;
ALTER table feature_type DROP KEY `name_class_idx`;
ALTER table feature_type ADD UNIQUE KEY `name_class_analysis_idx` (`name`,`class`, `analysis_id`);


#Need to alter enum here too
ALTER table feature_type MODIFY  `class` enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Regulatory Motif',  'Enhancer', 'Expression', 'Pseudo', 'Open Chromatin', 'Search Region', 'Association Locus', 'Segmentation State') default NULL;

analyze table feature_type;
optimize table feature_type;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_b.sql|feature_type.analysis_id');


