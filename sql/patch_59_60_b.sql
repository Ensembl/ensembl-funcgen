# patch_59_60_b.sql
#
# title: add pwm_feature
#
# description:
# Add pwm_feature table, and associated link table

DROP TABLE IF EXISTS `pwm_feature`;
CREATE TABLE `pwm_feature` (
  `pwm_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) default NULL,
  `score` double default NULL,
  `pwm_id` int(10) unsigned NOT NULL,	
  PRIMARY KEY  (`pwm_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `pwm_idx` (`pwm_id`),     
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;


DROP TABLE IF EXISTS `associated_pwm_feature`;
CREATE TABLE `associated_pwm_feature` (
   `annotated_feature_id` int(10) unsigned NOT NULL,
   `pwm_feature_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`annotated_feature_id`, `pwm_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;





-- pwm_feature to feature_set/cell_type may be 1 to many relationship
-- pwm_feature to annotated_feature may be many to 1 relationship
-- remove feature_set_id or change to cell_type_id?
-- cell_type style query will have to be done via query extension to feature_set using cell_type_id.




# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_58_59_b.sql|probe.description');


