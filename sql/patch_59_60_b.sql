# patch_59_60_b.sql
#
# title: associated_feature_type.table_name_id
#
# description:
# Change the feature_table/id to table_name/id

DROP TABLE IF EXISTS `new_associated_feature_type`;
CREATE TABLE `new_associated_feature_type` (
  `table_id` int(10) unsigned NOT NULL,
  `table_name` enum('annotated_feature','external_feature','regulatory_feature', 'feature_type') NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`table_id`,`table_name`,`feature_type_id`),
  KEY `feature_type_index` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



insert into new_associated_feature_type select feature_id, 'annotated_feature', feature_type_id from associated_feature_type where feature_table='annotated';
insert into new_associated_feature_type select feature_id, 'external_feature', feature_type_id from associated_feature_type where feature_table='external';
insert into new_associated_feature_type select feature_id, 'regulatory_feature', feature_type_id from associated_feature_type where feature_table='regulatory';


DROP table associated_feature_type;

CREATE TABLE `associated_feature_type` (
  `table_id` int(10) unsigned NOT NULL,
  `table_name` enum('annotated_feature','external_feature','regulatory_feature', 'feature_type') NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`table_id`,`table_name`,`feature_type_id`),
  KEY `feature_type_index` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


insert into associated_feature_type select * from new_associated_feature_type;

DROP table new_associated_feature_type;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_59_60_b.sql|associated_feature_type.table_name_id');


