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
@header patch_73_74_b.sql - input_set_subset_split
@desc   Patch to make input_subsets independant of input_sets
*/

INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'ChIP-Seq');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'DNase-Seq');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'FAIRE');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'RRBS');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'WGBS');

ALTER TABLE `input_set`    ADD `analysis_id`   smallint(5) unsigned NOT NULL;


UPDATE input_set inp, analysis a, feature_type ft SET inp.analysis_id = a.analysis_id 
  WHERE a.logic_name = 'ChIP-Seq' AND inp.feature_type_id = ft.feature_type_id 
  AND ft.class IN ('HISTONE', 'TRANSCRIPTION FACTOR', 'TRANSCRIPTION FACTOR COMPLEX', 'Polymerase');

UPDATE input_set inp, analysis a, analysis a1, result_set rs set inp.analysis_id=a.analysis_id
  WHERE a.logic_name='RRBS' AND a1.logic_name='RRBS_merged_filtered_10' 
  AND a1.analysis_id=rs.analysis_id AND rs.name=inp.name;
  
UPDATE input_set inp, analysis a, analysis a1, result_set rs set inp.analysis_id=a.analysis_id
  WHERE a.logic_name='WGBS' AND a1.logic_name='WGBS_merged' 
  AND a1.analysis_id=rs.analysis_id AND rs.name=inp.name;
  
  
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'DNase-Seq') WHERE `feature_type_id` IN ( SELECT `feature_type_id` FROM `feature_type` WHERE `name` = 'DNase1');
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'FAIRE')     WHERE `feature_type_id` IN ( SELECT `feature_type_id` FROM `feature_type` WHERE `name` = 'FAIRE');


UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.GM12878.comb11.concord4') WHERE `name` = 'Segmentation:GM12878';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.H1ESC.comb11.concord4')   WHERE `name` = 'Segmentation:H1ESC';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.HeLa-S3.comb11.concord4') WHERE `name` = 'Segmentation:HeLa-S3';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.HepG2.comb11.concord4')   WHERE `name` = 'Segmentation:HepG2';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.HUVEC.comb11.concord4')   WHERE `name` = 'Segmentation:HUVEC';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.K562.comb11.concord4')    WHERE `name` = 'Segmentation:K562';

DROP TABLE IF EXISTS `input_set_input_subset`;
CREATE TABLE `input_set_input_subset` (
    `input_set_id`    int(10) unsigned NOT NULL,
    `input_subset_id` int(10) unsigned NOT NULL,
    UNIQUE KEY `iset_subset_table_idx` (`input_subset_id`,`input_set_id`)
    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;

INSERT INTO input_set_input_subset ( input_set_id, input_subset_id ) SELECT input_set_id, input_subset_id FROM input_subset ;

ALTER TABLE `input_subset` ADD `cell_type_id`     int(10) unsigned DEFAULT NULL AFTER `input_subset_id`;
ALTER TABLE `input_subset` ADD `experiment_id`    int(10) unsigned NOT NULL AFTER `cell_type_id`;
ALTER TABLE `input_subset` ADD `feature_type_id`  int(10) unsigned NOT NULL AFTER `experiment_id`;

UPDATE input_set iset, input_subset iss SET iss.experiment_id = iset.experiment_id WHERE iss.input_set_id = iset.input_set_id;

ALTER TABLE input_subset DROP INDEX set_name_dx;
ALTER TABLE input_subset ADD UNIQUE name_exp_idx (name, experiment_id);

ALTER TABLE input_subset DROP COLUMN input_set_id;
ALTER TABLE input_set    DROP        format;
ALTER TABLE input_set    DROP        vendor;


UPDATE 
  input_subset isset, 
  input_set iset, 
  input_set_input_subset link 
SET 
  isset.feature_type_id = iset.feature_type_id, 
  isset.cell_type_id    = iset.cell_type_id 
WHERE 
  iset.input_set_id    = link.input_set_id     AND 
  link.input_subset_id = isset.input_subset_id;

-- ultimately input_set & input_set_input_subset will be dropped in favour of just using result_set
-- meaning input_set.analysis_id will move to input_subset. 
  
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_b.sql|input_set_subset_split');


