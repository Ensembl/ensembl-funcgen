--
-- Table structure for table `mage_xml`
--

DROP TABLE IF EXISTS `mage_xml`;
CREATE TABLE `mage_xml` (
   `mage_xml_id` int(10) unsigned NOT NULL auto_increment,
   `xml` text,
   PRIMARY KEY  (`mage_xml_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

alter table experiment add `mage_xml_id` int(10) unsigned default NULL;
alter table feature_set add `name` varchar(40) default NULL;
alter table result_set add `name` varchar(40) default NULL;
alter table data_set add `name` varchar(40) default NULL;


alter table experimental_chip change replicate `biological_replicate` varchar(40) default NULL;
alter table experimental_chip add `technical_replicate` varchar(40) default NULL;

--now add names to import result sets manually, exp_name_IMPORT
--coalesce Stunneburg chip rsets, making sure displayable is only set for the correct replicate
--populate experimental_chip replicate fields appropriately

alter table experiment change date `date` date default '0000-00-00';
update meta set meta_value=45 where meta_key='schema_version';
