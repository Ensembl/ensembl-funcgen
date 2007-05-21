--
-- Table structure for table `mage_xml`
--


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

-- add X and Y to result
alter table result add  `X` int(4) unsigned default NULL;
alter table result add  `Y` int(4) unsigned default NULL;


-- Need to update status for old ec and chans
-- select experimental_chip_id  from experimental_chip where experiment_id =12;

-- select e.name, rs.result_set_id from experiment e, result_set rs, experimental_chip ec, chip_channel cc where e.experiment_id=ec.experiment_id and ec.experiment_id>1 and ec.experiment_id <12 and ec.experimental_chip_id =cc.table_id and cc.table_name='experimental_chip' and cc.result_set_id=rs.result_set_id group by result_set_id;

create table `tmp_chip_channel`(
   `chip_channel_id` int(10) unsigned NOT NULL auto_increment,
   `result_set_id` int(10) unsigned default '0',
   `table_id` int(10) unsigned NOT NULL,
   `table_name` varchar(20) NOT NULL,
   PRIMARY KEY  (`result_set_id`, `chip_channel_id`),
   UNIQUE KEY `rset_table_idname_idx` (`result_set_id`, `table_id`, `table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


insert into tmp_chip_channel(select * from chip_channel);


---ctcf replicate chip_channel_id and name fix
---update chip_channel cc, tmp_chip_channel tcc set cc.chip_channel_id=tcc.chip_channel_id where tcc.table_name='experimental_chip' and cc.table_name='experimental_chip' and tcc.table_id=cc.table_id and tcc.result_set_id=17;
--update result_set set name=replace(name, 'SOM00H0', 'ctcf_ren');

drop table chip_channel;

rename table tmp_chip_channel to chip_channel;


-- rename biorep and techrep names as appropriate
--update result_set set name=replace(name, 'BIOREP', 'BR');
--update result_set set name=replace(name, 'techrep', 'TR');



alter table predicted_feature change display_label `display_label` varchar(60) default NULL;


---
-- Table structure for table `probe_design`
--

CREATE TABLE `probe_design` (
   `probe_id` int(10) unsigned NOT NULL default '0',
   `analysis_id` int(10) unsigned NOT NULL default '0',
   `score` double default NULL,	
   `coord_system_id` int(10) unsigned NOT NULL default '0',
    PRIMARY KEY  (`probe_id`, `analysis_id`, `coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



-- make feature/set names unique key
-- can't do result set yet due to duplicate names of channel and chip IMPORT sets
create unique index `name_idx` on feature_set (name);
create unique index `name_idx` on data_set (name);
