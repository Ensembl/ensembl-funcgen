/**
@header patch_73_74_b.sql - schema version
@desc   Update for new funcgen schema layout
*/
ALTER TABLE `status_name` MODIFY `name` varchar(60);

INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'ChIP-Seq');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'DNase-Seq');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'eQTL');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'FAIRE');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'PolIII');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'RRBS');
INSERT INTO `analysis` (`created`, `logic_name`) values (NOW(), 'WGBS');

ALTER TABLE `input_set`    ADD `analysis_id`   int(10) unsigned NOT NULL;


UPDATE
  `input_set`
SET
  `analysis_id` = (
    SELECT
      `analysis_id`
    FROM
      analysis
    WHERE
      `logic_name` = 'ChIP-Seq'
      )
WHERE
  `feature_type`_id` IN (
    SELECT
      `feature_type`_id`
    FROM
      `feature_type`
    WHERE
      `feature_type`.class IN ('HISTONE', 'TRANSCRIPTION FACTOR', 'TRANSCRIPTION FACTOR COMPLEX')
    )
  ;

UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'DNase-Seq') WHERE `feature_type`_id` IN ( SELECT `feature_type`_id` FROM `feature_type` WHERE `name` = 'DNase1';) ;
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'FAIRE')     WHERE `feature_type`_id` IN ( SELECT `feature_type`_id` FROM `feature_type` WHERE `name` = 'FAIRE';) ;
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'PolIII')    WHERE `feature_type`_id` IN ( SELECT `feature_type`_id` FROM `feature_type` WHERE `name` = 'PolIII';) ;

UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.GM12878.comb11.concord4') WHERE `name` = 'Segmentation:GM12878';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.H1ESC.comb11.concord4')   WHERE `name` = 'Segmentation:H1ESC';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.HeLa-S3.comb11.concord4') WHERE `name` = 'Segmentation:HeLa-S3';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.HepG2.comb11.concord4')   WHERE `name` = 'Segmentation:HepG2';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.HUVEC.comb11.concord4')   WHERE `name` = 'Segmentation:HUVEC';
UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'chromhmm.segway.K562.comb11.concord4')    WHERE `name` = 'Segmentation:K562';

UPDATE `input_set` SET `analysis_id` = ( SELECT `analysis_id` FROM `analysis` WHERE `logic_name` = 'eQTL')                                    WHERE `vendor` = 'EQTL';


DROP TABLE IF EXISTS `input_set_input_subset`;
CREATE TABLE `input_set_input_subset` (
    `input_set_id`    int(10) unsigned NOT NULL,
    `input_subset_id` int(10) unsigned NOT NULL,
    UNIQUE KEY `iset_subset_table_idx` (`input_subset_id`,`input_set_id`)
    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;

INSERT INTO input_set_input_subset ( input_set_id, input_subset_id ) SELECT input_set_id, input_subset_id FROM input_subset ;

ALTER TABLE `input_subset` ADD `cell_type_id`     int(10) unsigned DEFAULT NULL AFTER `input_subset_id`;
ALTER TABLE `input_subset` ADD `experiment_id`    int(10) unsigned DEFAULT NULL AFTER `cell_type_id`;
ALTER TABLE `input_subset` ADD `feature_type_id`  int(10) unsigned DEFAULT NULL AFTER `experiment_id`;

UPDATE input_set iset, input_subset iss SET iss.experiment_id = iset.experiment_id WHERE iss.input_set_id = iset.input_set_id;

ALTER TABLE input_subset DROP INDEX set_name_dx;
ALTER TABLE input_subset ADD UNIQUE name_exp_idx (name, experiment_id);

ALTER TABLE input_subset DROP COLUMN input_set_id;
ALTER TABLE `input_set`    DROP `format`;


