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

--- Data patch finished ---

-- alter xref table to capture coding and target info for external features

alter table xref change info_type `info_type` enum('PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED', 'CODING', 'TARGET') NOT NULL default 'PROJECTION';
update xref set info_type='TARGET' where info_type='DEPENDENT';


-- Add LOESS to status tables

INSERT into status_name values ('', 'LOESS');
-- Not really necessary, should add this to a store analysis script with a status flag


-- Add 'expression' class to feature_type
alter table feature_type modify class enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Overlap', 'Regulatory Motif', 'Region', 'Enhancer', 'Expression') default NULL;


-- Allow for longer names in experiment, experimental_chip, result_set
alter table experiment modify name varchar(100) default NULL;
alter table experimental_chip modify biological_replicate varchar(100) default NULL;
alter table experimental_chip modify technical_replicate varchar(100) default NULL;
alter table result_set modify name varchar(100) default NULL;
alter table experimental_subset modify name varchar(100) default NULL;



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
