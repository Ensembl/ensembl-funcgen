-- add table linking arrays to species

DROP TABLE IF EXISTS `array2organism`;

CREATE TABLE `array2organism` (
  `array2org_id` int(10) NOT NULL AUTO_INCREMENT,
  `array_id` int(27) NOT NULL,
  `organism` varchar(100) NOT NULL,
  PRIMARY KEY (`array2org_id`),
  KEY `array_id` (`array_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- add column to array table to specify whether array can be exported to
-- a species specific funcgen database

ALTER TABLE array ADD COLUMN  export2funcgen tinyint(4) NOT NULL DEFAULT '1';