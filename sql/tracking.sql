
DROP TABLE IF EXISTS `input_set_tracking`;
/* CREATE TABLE `input_set_tracking` (
  `input_set_id` int(10) unsigned NOT NULL,
  `release_version` varchar(45),
  `is_current` int(1),
  `status` enum('OK','REMOVE','REBUILD') default 'OK',
  PRIMARY KEY  (`input_set_id`),
  UNIQUE KEY `name_cell_type_feature_type_idx` (`experiment_name`,`species`,`cell_type`,`feature_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
*/

DROP TABLE IF EXISTS `result_set_stats`;

CREATE TABLE `result_set_stats` (
  `result_set_id` int(10) unsigned NOT NULL,
  `total_reads` int unsigned NOT NULL,
  `mapped_reads` int unsigned NOT NULL,
  `minimum_quality_reads` int unsigned NOT NULL,
  PRIMARY KEY  (`result_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- also add in feature_set_stats here!


DROP TABLE IF EXISTS `input_subset_tracking`;

CREATE TABLE `input_subset_tracking` (
  `input_subset_id` int(10) unsigned NOT NULL,
  `replicate` int unsigned default 1,
  `downloaded` datetime default NULL,
  `availability_date` datetime default NULL,
  `md5sum` varchar(45) default NULL,
  `local_url` text DEFAULT NULL,
  PRIMARY KEY  (`input_subset_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- add notes here as we have lost this from the input_set_tracking table?

-- Change this to current repository...
-- INSERT INTO meta (meta_key, meta_value) VALUES ('fastq_repository','/Your/species/specific/fastq/dir/');
-- INSERT INTO meta (meta_key, meta_value) VALUES ('current_coord_system','GRCh37');

-- Change this to something more appropriate tracking_only? (is_dev is ambiguous)
ALTER TABLE status_name ADD tracking_only TINYINT(1) NOT NULL default 0;

INSERT INTO status_name (name, tracking_only) values ('DOWNLOADED',              1);
INSERT INTO status_name (name, tracking_only) values ('IS_CURRENT',              1);
INSERT INTO status_name (name, tracking_only) values ('ADD_TO_REGULATORY_BUILD', 1);
INSERT INTO status_name (name, tracking_only) values ('IN_REGULATORY_BUILD',     1);
INSERT INTO status_name (name, tracking_only) values ('RELEASED',                1);
INSERT INTO status_name (name, tracking_only) values ('TO_BE_REVOKED',           1);
INSERT INTO status_name (name, tracking_only) values ('REVOKED',                 1);
INSERT INTO status_name (name, tracking_only) values ('TO_BE_REBUILD',           1);
INSERT INTO status_name (name, tracking_only) values ('REBUILT',                 1);
INSERT INTO status_name (name, tracking_only) values ('IN_RELEASE',              1);

-- what about status history?
-- We need to be able to track when a set was released, and whether it was pulled and then re-release
-- Separate table?


--ALTER TABLE `data_set_tracking`     DROP `is_current`;
DROP TABLE IF EXISTS `data_set_tracking`;
CREATE TABLE `data_set_tracking` (
  `data_set_id` int(10) unsigned NOT NULL,
  `release_version` varchar(45) DEFAULT NULL,
  `is_current` int(1) DEFAULT NULL,
  PRIMARY KEY (`data_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1; 


-- Do we really need data_set_tracking
-- This probably more accruately feature_set_tracking?
-- although we may also want to release stand alone result_sets





--ALTER TABLE `input_set_tracking`    DROP `status`;

