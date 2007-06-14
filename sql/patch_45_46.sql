-- add feature_type and cell_type columns result_set
alter table result_set add  `cell_type_id` int(10) unsigned default NULL;
alter table result_set add  `feature_type_id` int(10) unsigned default NULL;
drop index `name_analysis_idx` on result_set;
create unique index `unique_idx` on result_set (name, analysis_id, feature_type_id, cell_type_id);

-- now patch all cell_type and feature_type entries in result_set
--chip sets first
update result_set rs, chip_channel cc, experimental_chip ec set rs.feature_type_id=ec.feature_type_id where rs.result_set_id=cc.result_set_id and cc.table_name='experimental_chip' and cc.table_id=ec.experimental_chip_id;
update result_set rs, chip_channel cc, experimental_chip ec set rs.cell_type_id=ec.cell_type_id where rs.result_set_id=cc.result_set_id and cc.table_name='experimental_chip' and cc.table_id=ec.experimental_chip_id; 
-- now channel sets
update result_set rs, chip_channel cc, channel c, experimental_chip ec set rs.feature_type_id=ec.feature_type_id where rs.result_set_id=cc.result_set_id and cc.table_name='channel' and cc.table_id=c.channel_id and c.experimental_chip_id=ec.experimental_chip_id; 
update result_set rs, chip_channel cc, channel c, experimental_chip ec set rs.cell_type_id=ec.cell_type_id where rs.result_set_id=cc.result_set_id and cc.table_name='channel' and cc.table_id=c.channel_id and c.experimental_chip_id=ec.experimental_chip_id; 


-- change all ids to int(10), 
-- add UNIQUE KEY `name_idx` (name) on data_set
-- alter key on array_chip  UNIQUE KEY `array_design_idx` (`array_id`, `design_id`)
-- recreate predicted_feature as annotated_feature

-- add enum on channel type TOTAL, EXPERIMENTAL & psuedo channels

-- patch feature set for ctcf, re-add new nessie analysis 14



-- add core regulatory tables as supporting_feature tables
-- regulatory_factor_coding is empty and unused?


-- supporting feature table or import directly into annotated feature
-- where do we draw the line between what goes in supporting rather than annotated?
-- what do we do about the overloading of the feature_type table?
-- supporting features must have multiple insatnces of feature_type with unique ids
-- i.e. high volume e.g. individual miRNAs

CREATE TABLE `supporting_feature` (
 `supporting_feature_id` int(10) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `analysis_id` smallint(5) unsigned NOT NULL default '0',
  `regulatory_factor_id` int(10) unsigned default NULL,
  `coord_system_id` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`supporting_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`analysis_id`,`seq_region_start`),
  KEY `seq_region_idx_2` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- imported core regulatory_feature/factor/search_region here
--insert into supporting_feature(select *, 1 from regulatory_feature);

--strip off ENST prefixes for names??


-- add feature_type_id column to replace type
alter table regulatory_factor add feature_type_id int(10) unsigned NOT NULL default '0';
insert into feature_type(name, class, description) values('miRNA Target', 'RNA', 'miRNA target motif');
insert into feature_type(name, class, description) values('Transcription Factor', 'TRANSCRIPTION FACTOR', 'Transcription factor motif');
insert into feature_type(name, class, description) values('Transcription Factor Complex', 'TRANSCRIPTION FACTOR', 'Transcription complex factor motif');

update regulatory_factor rf, feature_type ft set rf.feature_type_id=ft.feature_type_id where rf.type='miRNA_target' and ft.name='miRNA Target';

--don't need to update other as they are all NULL

-- remove old type column

alter table regulatory_factor drop type;


--we need to split the regulatory_search_region table to extract the xrefs


-- add regulatory factor types to feature_type
