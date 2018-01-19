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
@header patch_68_69_b.sql - DNAMethylationFeature support
@desc   add a feature_class field to result_set
*/


-- Move currently misplaced fields from experiment table

ALTER TABLE result_set ADD `feature_class` enum('result','dna_methylation') DEFAULT NULL;
ALTER TABLE result_set drop index `unique_idx`;
ALTER TABLE result_set ADD UNIQUE INDEX `unique_idx` (`name`,`analysis_id`,`feature_type_id`,`cell_type_id`, `feature_class`);


ALTER TABLE feature_type MODIFY class enum('Insulator','DNA','Regulatory Feature','Histone','RNA','Polymerase','Transcription Factor','Transcription Factor Complex','Regulatory Motif','Enhancer','Expression','Pseudo','Open Chromatin','Search Region','Association Locus','Segmentation State', 'DNA Modification') DEFAULT NULL;

-- actully want to remove input_set feature class as input_set can support more than one feature class

ALTER TABLE input_set MODIFY type  enum('annotated','result','segmentation','dna_methylation') DEFAULT NULL;

UPDATE result_set set feature_class = 'result';

OPTIMIZE TABLE result_set;
ANALYZE TABLE result_set;

OPTIMIZE TABLE feature_type;
ANALYZE TABLE feature_type;


OPTIMIZE TABLE input_set;
ANALYZE TABLE input_set;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_68_69_b.sql|DNAMethylationFeature support');


