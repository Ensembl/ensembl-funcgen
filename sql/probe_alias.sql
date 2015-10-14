/* add table for probe alias tracking - /*
/* column default_name =1 indicates this entry is the default name for the probe*/

DROP TABLE IF EXISTS `probe_alias`;
CREATE TABLE `probe_alias` (
       `probe_alias_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
       `probe_id` int(10) NOT NULL,
       `alias` varchar(100) NOT NULL,
       `default_name` BOOL DEFAULT 0 ,
       PRIMARY KEY (`probe_alias_id`),
       UNIQUE KEY `probe_sha1_idx` (`probe_id`, `alias`),
       KEY `alias_idx` (`alias`),
       KEY `alias_probe_id` (`probe_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

