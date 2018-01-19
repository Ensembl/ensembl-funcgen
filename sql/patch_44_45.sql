-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

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

--create table `tmp_chip_channel`(
 --  `chip_channel_id` int(10) unsigned NOT NULL auto_increment,
 --  `result_set_id` int(10) unsigned default '0',
 --  `table_id` int(10) unsigned NOT NULL,
 --  `table_name` varchar(20) NOT NULL,
 --  PRIMARY KEY  (`result_set_id`, `chip_channel_id`),
 --  UNIQUE KEY `rset_table_idname_idx` (`result_set_id`, `table_id`, `table_name`)
--) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--this primary key is not increenting the cc_id when we are using a different table_name
create table `tmp_chip_channel`(
   `chip_channel_id` int(10) unsigned NOT NULL auto_increment,
   `result_set_id` int(10) unsigned default '0',
   `table_id` int(10) unsigned NOT NULL,
   `table_name` varchar(20) NOT NULL,
   PRIMARY KEY  (`chip_channel_id`, `result_set_id`),
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
-- can't do this until we split the name of into a separate table

--create unique index `name_idx` on data_set (name);


--insert into cell_type values(NULL, 'CD4', NULL, 'Human CD4 T-Cells');
--update experimental_chip ec, cell_type ct set ec.cell_type_id=ct.cell_type_id where ct.name='CD4' and ec.unique_id='CD4_parzen_02';
--update experimental_chip ec, cell_type ct set ec.cell_type_id=ct.cell_type_id where ct.name='GM06996' and ec.unique_id='GM06990_parzen_0115';

--insert into feature_type values('', 'DNase', 'DNA', 'DNase Hypersensitive Site');
--update experimental_chip ec, feature_type ft set ec.feature_type_id=ft.feature_type_id where ft.name='DNase' and ec.unique_id='CD4_parzen_02';
--update experimental_chip ec, feature_type ft set ec.feature_type_id=ft.feature_type_id where ft.name='DNase' and ec.unique_id='GM06990_parzen_0115';



-- update  result_set rs, data_set ds set ds.name=rs.name where rs.result_set_id=ds.result_set_id;
-- update  feature_set fs, data_set ds set fs.name=ds.name where fs.feature_set_id=ds.feature_set_id;
-- insert into feature_type values ('', 'CTCF', 'INSULATOR', 'CCCCTC-binding factor');
-- update experimental_chip ec, experiment e, feature_type ft  set ec.feature_type_id=ft.feature_type_id where e.name='ctcf_ren' and e.experiment_id=ec.experiment_id and ft.name='CTCF';
-- insert into cell_type values ('', 'IMR90', '', 'Human Fetal Lung Fibroblast');
-- update cell_type set display_label =NULL where display_label='';
-- update experimental_chip ec, experiment e, cell_type ct  set ec.cell_type_id=ct.cell_type_id where e.name='ctcf_ren' and e.experiment_id=ec.experiment_id and ct.name='IMR90';



-- insert into cell_type values(NULL, 'HL-60', NULL, 'Hman promyelotic Leukemia Cells');
--  update feature_set set cell_type_id=6 where name like "Wiggle%";

drop index analysis_idx on result_set;
create unique index `name_analysis_idx` on result_set (name, analysis_id);


delete s from status s, status_name sn where sn.name='DISPLAYABLE' and sn.status_name_id=s.status_name_id;
insert into status(table_id, table_name, status_name_id) values(24, 'data_set', 2);
insert into status(table_id, table_name, status_name_id) values(27, 'feature_set', 2);

-- patch feature_set table to allow for longer names
ALTER TABLE feature_set CHANGE name `name` varchar(250) DEFAULT NULL;



--set on RegulatoryFeatures to displayble
delete s from status s, status_name sn where sn.name='DISPLAYABLE' and sn.status_name_id=s.status_name_id;
insert into status(table_id, table_name, status_name_id) values(24, 'data_set', 2);
insert into status(table_id, table_name, status_name_id) values(27, 'feature_set', 2);

--feature class update to facilitate get data set methods
update feature_type set class='REGULATORY FEATURE' where feature_type_id =19;

--feature type correction
update feature_type set name=replace(name, 'H3K20', 'H4K20');
update feature_set set name=replace(name, 'H3K20', 'H4K20');
update feature_type set class='HISTONE' where feature_type_id in(11,12,13,14); 
update feature_type set class='OVERLAP' where feature_type_id in(16, 17,18); 
update data_set set name ='Wiggle_H4K20me3' where name='H3K20me3_wiggle';



--CTCF displayable update
insert into status(table_id, table_name, status_name_id) select experimental_chip_id, 'experimental_chip', 2 from experimental_chip where experiment_id =12;
insert into status values(22, 'feature_set', 2);
insert into status values(22, 'data_set', 2);
insert into status values(18, 'result_set', 2);

--remove CTCF spurious data set and revert rset name back to original
delete from data_set where data_set_id=23;
update result_set set name = 'ctcf_ren_BR1' where name ='ctcf_ren_BR1_TR2';

-- correct ctcf analysis id
update feature_set set analysis_id =7 where analysis_id =14;



--- unique index on exp name
create unique index `name_idx` on experiment (name);


