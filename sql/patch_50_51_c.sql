
--- Dev stuff





DROP TABLE IF EXISTS `associated_feature_type`;
CREATE TABLE `associated_feature_type` (
   `feature_id` int(10) unsigned NOT NULL,
   `feature_table` enum('annotated', 'external', 'regulatory') default NULL,
   `feature_type_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`feature_id`, `feature_table`, `feature_type_id`),
   KEY `feature_type_index` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_c.sql|associated_feature_type');
