select "This patch should be worked through manually as there are data specific and species specific patches";
exit;



-- Human Data patch
-- Fix absent feature_type info
update feature_type set class ='Histone', description='Histone 4 Lysine 20 Tri-Methylation' where name='H4K20me3';
update feature_type set description='Histone 3 Lysine 27 Tri-Methylation' where name='H3K27me3';
update feature_type set description='Histone 3 Lysine 36 Tri-Methylation' where name='H3K36me3';
update feature_type set description='Histone 3 Lysine 9 Tri-Methylation' where name='H3K9me3';
update feature_type set description='Histone 3 Lysine 79 Tri-Methylation' where name='H3K79me3';
--tidy up array description
update array set description='2005-05-10_HG17Tiling_Set. Whole human genome (hg17 from UCSC) tiled at 100 bp spacing. Repeat masked. 148 synthesis cycle limit. Contains random probes with GC content between 15-35%.' where name='2005-05-10_HG17Tiling_Set';
--remove old ChiP-Seq pseudo array_chip
-- Data patch finished




--change status_name to remove IMPORTED from analysis states
--Need to change this in the API after branch-48!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
update status_name set name=replace(name, 'IMPORTED_', '');


--migrate supporting_set_type to supporting_set table, to allow multiple supporting_set type for the same data_set
-- also remove external from enum as this is a feature_set
--NEED TO CHANGE IN API and SQL AFTER v48 branch!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
alter table supporting_set add `type` enum('result','feature','experimental') default NULL;
alter table supporting_set add key `type_idx` (type);
update supporting_set ss, data_set ds set ss.type=ds.supporting_set_type where ss.data_set_id=ds.data_set_id;
alter table data_set drop column supporting_set_type;


--change egroup to experimental_group
--NEED TO CHANGE THIS IN THE API and the SQL AFTER THE V48 BRANCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

--change DAS DISPLAYABLE to DAS_DISPLAYABLE?


--Need to patch supporting_set, migrate type from data_set to allow multiple types of supporting set

--change channel to enum CONTROL/EXPERIMETNAL can we not NULL this also? this will default ot 1st no?


--compare mysql median to perl median for ResultFeature query

--need to vhange average row length on this and on regulatory feature!!



-- tidy up overlap feature_sets and create data_sets for them?
-- reduce size of name field in feature_set


-- should change max rows to 17000000 for reg feats to reflect max stable ids
-- also consider average row?
-- change small table primary key ids to medium int?



--add key on ec cell_type_id?
-- add enum on channel type TOTAL, EXPERIMENTAL & DUMMY? channels

-- add core regulatory tables as other_feature tables
-- regulatory_factor_coding is empty and unused? migrate to xrefs


--enum feature_type class default NULL
