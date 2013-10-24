
DROP TABLE IF EXISTS `input_set_tracking`;
CREATE TABLE `input_set_tracking` (
  `input_set_id` int(10) unsigned NOT NULL,
  `release_version` varchar(45),
  `is_current` int(1),
  `status` enum('OK','REMOVE','REBUILD') default 'OK',
  PRIMARY KEY  (`input_set_id`),
  UNIQUE KEY `name_cell_type_feature_type_idx` (`experiment_name`,`species`,`cell_type`,`feature_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `result_set_stats`;

CREATE TABLE `result_set_stats` (
  `result_set_id` int(10) unsigned NOT NULL,
  `total_reads` int unsigned NOT NULL,
  `mapped_reads` int unsigned NOT NULL,
  `minimum_quality_reads` int unsigned NOT NULL,
  PRIMARY KEY  (`result_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `input_subset_tracking`;

CREATE TABLE `input_subset_tracking` (
  `input_subset_id` int(10) unsigned NOT NULL,
  `replicate` int unsigned default 1,
  `downloaded` datetime default NULL,
  `availability_date` datetime default NULL,
  `md5sum` varchar(45) default NULL,
  PRIMARY KEY  (`input_subset_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Change this to current repository...
INSERT INTO meta (meta_key, meta_value) VALUES ('repository','/lustre/scratch109/ensembl/funcgen/tj1/pipeline/data_home/fastq');

ALTER TABLE status_name ADD is_dev TINYINT(1) NOT NULL default 0;

INSERT INTO status_name (name) values ('DOWNLOADED', 0);
INSERT INTO status_name (name) values ('IS_CURRENT', 0);
INSERT INTO status_name (name) values ('ADD_TO_REGULATORY_BUILD', 0);
INSERT INTO status_name (name) values ('IN_REGULATORY_BUILD', 0);
INSERT INTO status_name (name) values ('RELEASED', 0);
INSERT INTO status_name (name) values ('TO_BE_REVOKED', 0);
INSERT INTO status_name (name) values ('REVOKED', 0);
INSERT INTO status_name (name) values ('TO_BE_REBUILD', 0);
INSERT INTO status_name (name) values ('REBUILT', 0);
INSERT INTO status_name (name) values ('IN_RELEASE', 0);


ALTER TABLE `data_set_tracking`     DROP `is_current`;
ALTER TABLE `input_set_tracking`    DROP `status`;
ALTER TABLE `input_subset_tracking` ADD `local_url` text DEFAULT NULL;

INSERT INTO meta (meta_key, meta_value) VALUES ('current_coord_system','GRCh37');
