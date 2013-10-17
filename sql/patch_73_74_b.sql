/**
@header patch_73_74_b.sql - schema version
@desc   Update for new funcgen schema layout
*/
ALTER TABLE status_name MODIFY name varchar(60);

INSERT INTO status_name (name) values ('IS_CURRENT');
INSERT INTO status_name (name) values ('DOWNLOADED');
INSERT INTO status_name (name) values ('IS_CONTROL');
INSERT INTO status_name (name) values ('ADD_TO_REGULATORY_BUILD');
INSERT INTO status_name (name) values ('IN_REGULATORY_BUILD');
INSERT INTO status_name (name) values ('RELEASED');
INSERT INTO status_name (name) values ('REVOKE');
INSERT INTO status_name (name) values ('REVOKED');
INSERT INTO status_name (name) values ('REBUILD');
INSERT INTO status_name (name) values ('REBUILT');

INSERT INTO analysis (created, logic_name) values (NOW(), 'ChIP-Seq');
INSERT INTO analysis (created, logic_name) values (NOW(), 'eQTL');
INSERT INTO analysis (created, logic_name) values (NOW(), 'RRBS');
INSERT INTO analysis (created, logic_name) values (NOW(), 'WGBS');
INSERT INTO analysis (created, logic_name) values (NOW(), 'FAIRE');

ALTER TABLE `input_set`    ADD `analysis_id`   int(10) unsigned NOT NULL;

UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'chromhmm.segway.GM12878.comb11.concord4' and iset.name = 'Segmentation:GM12878';
UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'chromhmm.segway.H1ESC.comb11.concord4'   and iset.name = 'Segmentation:H1ESC';
UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'chromhmm.segway.HeLa-S3.comb11.concord4' and iset.name = 'Segmentation:HeLa-S3';
UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'chromhmm.segway.HepG2.comb11.concord4'   and iset.name = 'Segmentation:HepG2';
UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'chromhmm.segway.HUVEC.comb11.concord4'   and iset.name = 'Segmentation:HUVEC';
UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'chromhmm.segway.K562.comb11.concord4'    and iset.name = 'Segmentation:K562';
UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'eQTL'                                    and iset.vendor = 'EQTL';
UPDATE input_set iset, analysis a SET iset.analysis_id = a.analysis_id WHERE a.logic_name = 'ChIP-Seq'                                and iset.analysis_id is NULL;

ALTER TABLE `status` ADD `timestamp` datetime DEFAULT NULL;
ALTER TABLE `status` ADD `release` int(4)     DEFAULT NULL;


ALTER TABLE `data_set_tracking`     DROP `is_current`;
ALTER TABLE `input_set_tracking`    DROP `status`;
ALTER TABLE `input_subset_tracking` ADD `local_url` text DEFAULT NULL;

DROP TABLE IF EXISTS `input_set_input_subset`;
CREATE TABLE `input_set_input_subset` (
    `input_set_id`    int(10) unsigned NOT NULL,
    `input_subset_id` int(10) unsigned NOT NULL,
    UNIQUE KEY `iset_subset_table_idx` (`input_subset_id`,`input_set_id`)
    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;

INSERT INTO input_set_input_subset ( input_set_id, input_subset_id ) SELECT input_set_id, input_subset_id FROM input_subset ;
ALTER TABLE `input_subset` ADD `experiment_id` int(10) unsigned DEFAULT NULL AFTER `input_subset_id`;
UPDATE input_set iset, input_subset iss SET iss.experiment_id = iset.experiment_id WHERE iss.input_set_id = iset.input_set_id;

ALTER TABLE input_subset DROP INDEX set_name_dx;
ALTER TABLE input_subset ADD UNIQUE name_exp_idx (name, experiment_id);

ALTER TABLE input_subset DROP COLUMN input_set_id;
ALTER TABLE `input_set`    DROP `format`;


