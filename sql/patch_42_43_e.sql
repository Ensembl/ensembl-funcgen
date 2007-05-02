-- create chip_channel table and entries from result
-- add chip_channel_id to result
-- drop table_id/name index on result and create new one on chip_channel_id
-- create result_set and entries from experimental_chip and update chip_channel result_set_ids

-- create chip_channel and entries from result


select "There are manual steps to this patch dependent on the species DB you are working with";
exit;

CREATE TABLE `chip_channel` (
   `chip_channel_id` int(10) unsigned NOT NULL auto_increment,
   `result_set_id` int(10) unsigned default '0',
   `table_id` int(10) unsigned  NULL,
   `table_name` varchar(20) NOT NULL,
   PRIMARY KEY  (`chip_channel_id`),
   UNIQUE KEY `rset_table_idname_idx` (`result_set_id`, `table_id`, `table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into chip_channel (table_id, table_name) select distinct(table_id), table_name from result;
alter table result add  `chip_channel_id` int(10) unsigned default NULL;
update result r, chip_channel cc set r.chip_channel_id=cc.chip_channel_id where cc.table_name=r.table_name and cc.table_id=r.table_id;
create INDEX chip_channel_idx on result (chip_channel_id);
drop INDEX table_name_id_idx on result;
alter table result drop table_id;
alter table result drop table_name;
alter table result change  chip_channel_id chip_channel_id int(10) unsigned NOT NULL;



-- manually add replicate info for Nimblegen stuff
-- for sanger rstuff, just add all ec's from same epxeriment to same result set


DROP TABLE IF EXISTS `result_set`;
CREATE TABLE `result_set` (
   `result_set_id` int(10) unsigned NOT NULL auto_increment,
   `analysis_id` int(10) unsigned default NULL,
   PRIMARY KEY  (`result_set_id`),
   KEY  `analysis_idx` (`analysis_id`) 
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--manually add result_sets
--insert into result_set(analysis_id) select r.analysis_id from  result r, experimental_chip ec, chip_channel cc where ec.replicate='UNKNOWN' and ec.experimental_chip_id=cc.table_id and cc.table_name='experimental_chip' and cc.chip_channel_id=r.chip_channel_id group by ec.experiment_id;
--update chip_channel cc,  experimental_chip ec set cc.result_set_id=(ec.experiment_id - 1) where ec.replicate='UNKNOWN' and ec.experimental_chip_id=cc.table_id and cc.table_name='experimental_chip';
--created 3 more result sets corresponding to channel analysis and two nimblegen chip analyses
--update chip_channel cc, experimental_chip ec set cc.result_set_id=12 where cc.table_name='experimental_chip' and cc.table_id=ec.experimental_chip_id and ec.replicate=1;


alter table result drop analysis_id;
