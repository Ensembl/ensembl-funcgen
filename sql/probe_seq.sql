--
-- Table structure for table `probe_seq`
--

DROP TABLE IF EXISTS `probe_seq`;

CREATE TABLE `probe_seq` (
  `probe_seq_id` int(10) NOT NULL AUTO_INCREMENT,
  `probe_sha1` char(40) NOT NULL,
  `probe_dna` text NOT NULL,
  `has_been_mapped` BOOL DEFAULT 0,
  PRIMARY KEY (`probe_seq_id`),
  UNIQUE KEY `probe_sha1_idx` (`probe_sha1`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

ALTER TABLE probe ADD COLUMN probe_seq_id int(10) DEFAULT NULL;

ALTER TABLE probe ADD KEY probe_seq_id (`probe_seq_id`);


DROP TABLE IF EXISTS `probe_alias`;
CREATE TABLE `probe_alias` (
       `probe_alias_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
       `probe_id` int(10) NOT NULL,
       `alias` varchar(100) NOT NULL,
       `default_name` tinyint(1) NOT NULL DEFAULT 0 ,
       PRIMARY KEY (`probe_alias_id`),
       KEY `alias_idx` (`alias`),
       KEY `alias_probe_id` (`probe_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
