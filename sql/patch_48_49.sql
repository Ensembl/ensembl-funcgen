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

select "This patch should be worked through manually as there are data specific and species specific patches";
exit;

--- MAKE SURE YOU EDIT efg.sql AFTER ADDING A PATCH ! ---

--- Human Data patch ---

-- Fix absent feature_type info
update feature_type set class ='Histone', description='Histone 4 Lysine 20 Tri-Methylation' where name='H4K20me3';
update feature_type set description='Histone 3 Lysine 27 Tri-Methylation' where name='H3K27me3';
update feature_type set description='Histone 3 Lysine 36 Tri-Methylation' where name='H3K36me3';
update feature_type set description='Histone 3 Lysine 9 Tri-Methylation' where name='H3K9me3';
update feature_type set description='Histone 3 Lysine 79 Tri-Methylation' where name='H3K79me3';

--tidy up array description
update array set description='2005-05-10_HG17Tiling_Set. Whole human genome (hg17 from UCSC) tiled at 100 bp spacing. Repeat masked. 148 synthesis cycle limit. Contains random probes with GC content between 15-35%.' where name='2005-05-10_HG17Tiling_Set';

--Remove old Human ChIP-Seq ProbeFeature/Probe/ArrayChip
--Was missed due to truncation of cell type name in design_id
delete pf from probe_feature pf, probe p, array_chip ac where ac.design_id='GM06990_DNASE:GM0699' and p.array_chip_id=ac.array_chip_id and pf.probe_id=p.probe_id;
delete p from probe p, array_chip ac where ac.design_id='GM06990_DNASE:GM0699' and p.array_chip_id=ac.array_chip_id;
delete from array_chip where design_id='GM06990_DNASE:GM0699';

--update feature_type table, add missing class and description fields
update feature_type set class='Histone' where name='H4K20me3';
update feature_type set description='Histone 4 Lysine 20 Tri-Methylation' where name='H4K20me3';
update feature_type set description='Histone 3 Lysine 27 Tri-Methylation' where name='H3K27me3';
update feature_type set description='Histone 3 Lysine 36 Tri-Methylation' where name='H3K36me3';
update feature_type set description='Histone 3 Lysine 79 Tri-Methylation' where name='H3K79me3';
update feature_type set description='Histone 3 Lysine 9 Tri-Methylation' where name='H3K9me3';

--- Human Data patch finished ---

-- alter xref table to capture coding and target info for external features

alter table xref change info_type `info_type` enum('PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED', 'CODING', 'TARGET') NOT NULL;
update xref set info_type='TARGET' where info_type='DEPENDENT';
alter table object_xref modify `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature', 'FeatureType') NOT NULL;


-- Add LOESS to status tables

INSERT into status_name values ('', 'LOESS');
-- Not really necessary, should add this to a store analysis script with a status flag


-- Add 'expression' class to feature_type
alter table feature_type modify class enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Overlap', 'Regulatory Motif', 'Region', 'Enhancer', 'Expression') default NULL;


-- Allow for longer names in experiment, experimental_chip, *_set, and analysis
alter table experiment modify name varchar(100) default NULL;
alter table experimental_chip modify biological_replicate varchar(100) default NULL;
alter table experimental_chip modify technical_replicate varchar(100) default NULL;
alter table result_set modify name varchar(100) default NULL;
alter table data_set modify name varchar(100) default NULL;
alter table feature_set modify name varchar(100) default NULL;
alter table experimental_subset modify name varchar(100) default NULL;
alter table analysis modify logic_name varchar(100) default NULL;

-- Some minimal filed definition tidy up
-- Most of these are due to differences in table creations between MySQL 4 > 5
-- NOT NULLs where extended to default '' or 0 dependent on variable.
-- enums no longer extended to explicitly state default logic
alter table array_chip modify `array_id` int(10) unsigned NOT NULL;
alter table cell_type modify `name` varchar(120) NOT NULL;
alter table channel modify `experimental_chip_id` int(10) unsigned default NULL;
alter table channel drop description;
alter table chip_channel modify `table_id` int(10) unsigned NOT NULL;
alter table chip_channel modify `table_name` varchar(20) NOT NULL;
alter table chip_channel modify `result_set_id` int(10) unsigned NOT NULL;
alter table coord_system modify `core_coord_system_id` int(10) NOT NULL;
alter table experimental_chip modify `experiment_id` int(10) unsigned default NULL;
alter table experimental_chip modify `array_chip_id` int(10) unsigned default NULL;
alter table experimental_set modify  `name` varchar(40) not NULL;
alter table experimental_subset modify `experimental_set_id` int(10) unsigned NOT NULL;
alter table external_db modify `db_name` varchar(28) NOT NULL;
alter table external_db modify status ENUM('KNOWNXREF','KNOWN','XREF','PRED','ORTH', 'PSEUDO') NOT NULL;
alter table external_db modify `priority` int(11) NOT NULL;
alter table external_synonym modify `xref_id` int(10) unsigned NOT NULL;
alter table external_synonym modify `synonym` varchar(40) NOT NULL;
alter table feature_set modify `feature_type_id` int(10) unsigned NOT NULL;
alter table feature_type modify  `name` varchar(40) NOT NULL;
alter table go_xref modify linkage_type ENUM('IC', 'IDA', 'IEA', 'IEP', 'IGI', 'IMP', 'IPI', 'ISS', 'NAS', 'ND', 'TAS', 'NR', 'RCA') NOT NULL;

alter table identity_xref modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table identity_xref modify `object_xref_id` int(10) unsigned NOT NULL;
alter table object_xref modify `object_xref_id` int(10) unsigned NOT NULL auto_increment;
alter table object_xref modify  `ensembl_id` int(10) unsigned NOT NULL;
alter table object_xref modify  `xref_id` int(10) unsigned NOT NULL;
alter table result modify `chip_channel_id` int(10) unsigned NOT NULL;

alter table status modify `status_name_id` int(10) NOT NULL;
alter table status_name modify   `name` varchar(20) default NULL;
alter table xref modify `dbprimary_acc` varchar(40) NOT NULL;
alter table xref modify   `display_label` varchar(128) NOT NULL;
alter table xref modify  `external_db_id` smallint(5) unsigned NOT NULL;

alter table annotated_feature modify `seq_region_id` int(10) unsigned NOT NULL;
alter table annotated_feature modify `seq_region_start` int(10) unsigned NOT NULL;
alter table annotated_feature modify `seq_region_end` int(10) unsigned NOT NULL;
alter table annotated_feature modify   `seq_region_strand` tinyint(1) NOT NULL;
alter table annotated_feature modify `feature_set_id` int(10) unsigned NOT NULL;

alter table external_feature modify  `seq_region_id` int(10) unsigned NOT NULL;
alter table external_feature modify   `seq_region_start` int(10) unsigned NOT NULL;
alter table external_feature modify   `seq_region_end` int(10) unsigned NOT NULL;
alter table external_feature modify   `seq_region_strand` tinyint(1) NOT NULL;
alter table external_feature modify  `feature_set_id` int(10) unsigned NOT NULL;

alter table regulatory_feature modify  `seq_region_id` int(10) unsigned NOT NULL;
alter table regulatory_feature modify   `seq_region_start` int(10) unsigned NOT NULL;
alter table regulatory_feature modify   `seq_region_end` int(10) unsigned NOT NULL;
alter table regulatory_feature modify   `seq_region_strand` tinyint(1) NOT NULL;

alter table regulatory_attribute modify  `regulatory_feature_id` int(10) unsigned NOT NULL;
alter table regulatory_attribute modify  `attribute_feature_id` int(10) unsigned NOT NULL;

alter table experimental_chip modify `unique_id` varchar(20) NOT NULL;
alter table analysis_description modify `analysis_id` int(10) unsigned NOT NULL;

alter table meta modify `meta_key` varchar(40) NOT NULL;
alter table meta modify `meta_value` varchar(255) NOT NULL;
alter table meta_coord modify `table_name` varchar(40) NOT NULL;
alter table meta_coord modify `coord_system_id` int(10) NOT NULL;
alter table coord_system modify  `name` varchar(40) NOT NULL;
alter table coord_system modify   `rank` int(11) NOT NULL;
alter table coord_system modify `coord_system_id` int(10) NOT NULL auto_increment;

alter table seq_region modify `name` varchar(40) NOT NULL;
alter table seq_region modify `coord_system_id` int(10) unsigned NOT NULL;
alter table seq_region modify `core_seq_region_id` int(10) unsigned NOT NULL;

alter table probe_feature modify  `seq_region_id` int(10) unsigned NOT NULL;
alter table probe_feature modify   `seq_region_start` int(10) NOT NULL;
alter table probe_feature modify   `seq_region_end` int(10) NOT NULL;
alter table probe_feature modify   `seq_region_strand` tinyint(4) NOT NULL;
alter table probe_feature modify   `probe_id` int(10) unsigned NOT NULL;
alter table probe_feature modify   `analysis_id` int(10) unsigned NOT NULL;
alter table probe_feature modify   `mismatches` tinyint(4) NOT NULL;

alter table probe_set modify `name` varchar(20) NOT NULL;
alter table probe_set modify `size` smallint(6) unsigned NOT NULL;

alter table probe modify  `name` varchar(40) NOT NULL;
alter table probe modify  `length` smallint(6) unsigned NOT NULL;
alter table probe modify  `array_chip_id` int(10) unsigned NOT NULL;

alter table probe_design modify `probe_id` int(10) unsigned NOT NULL;
alter table probe_design modify `analysis_id` int(10) unsigned NOT NULL;
alter table probe_design modify `coord_system_id` int(10) unsigned NOT NULL;

alter table supporting_set modify `data_set_id` int(10) unsigned NOT NULL;
alter table supporting_set modify `supporting_set_id` int(10) unsigned NOT NULL;


CREATE TABLE `experimental_design` (
   `design_type_id` int(10) unsigned NOT NULL auto_increment,
   `table_name` varchar(40) default NULL,
   `table_id` int(10) unsigned default NULL,	
   PRIMARY KEY  (`design_type_id`, `table_name`, `table_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- likley no data in here anyway
insert into experimental_design select * from experimental_design_type;
DROP TABLE IF EXISTS `experimental_design_type`;










--change status_name to remove IMPORTED from analysis states
--NEED TO TEST
update status_name set name=replace(name, 'IMPORTED_', '');


--migrate supporting_set_type to supporting_set table, to allow multiple supporting_set type for the same data_set
-- also remove external from enum as this is a feature_set
--NEED TO TEST!
alter table supporting_set add `type` enum('result','feature','experimental') default NULL;
alter table supporting_set add key `type_idx` (type);
update supporting_set ss, data_set ds set ss.type=ds.supporting_set_type where ss.data_set_id=ds.data_set_id;
alter table data_set drop column supporting_set_type;
--WE NEED TO DROP THE KEY TOO!


--change egroup to experimental_group
CREATE TABLE `experimental_group` (
  `experimental_group_id` smallint(6) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL default '',
  `location` varchar(120) default NULL,
  `contact` varchar(40) default NULL,
  PRIMARY KEY  (`experimental_group_id`),
  UNIQUE KEY name_idx(`name`)
) ENGINE=MyISAM AUTO_INCREMENT=2 DEFAULT CHARSET=latin1;

insert into experimental_group select * from egroup;
drop table egroup;
alter table experiment change egroup_id experimental_group_id smallint(6) unsigned default NULL;
alter table experiment drop key egroup_idx;
alter table experiment add KEY `experimental_group_idx` (`experimental_group_id`);


--change DAS DISPLAYABLE to DAS_DISPLAYABLE
update status_name set name='DAS_DISPLAYABLE' where name ='DAS DISPLAYABLE';


-- Create exp_chip IMPORT set for Stunnenburg set and create toplebel and BR TR sets?
-- Not essential

--drop description from channel?
--Can we clean DUMMY channel entries by changing the sql to a left join?
--what was the problem here, ResultSet adadptor?



-- add enum on channel type TOTAL, EXPERIMENTAL & DUMMY? channels
--change channel to enum CONTROL/EXPERIMETNAL can we not NULL this also? this will default ot 1st no?


--compare mysql median to perl median for ResultFeature query

--need to change average row length on this and on regulatory feature!!



-- tidy up overlap feature_sets and create data_sets for them?
-- reduce size of name field in feature_set


-- should change max rows to 17000000 for reg feats to reflect max stable ids? Is this sufficient for at least two regulatory builds?
-- also consider average row length?
-- change small table primary key ids to medium int?



--add key on ec cell_type_id?
--enum feature_type class default NULL
