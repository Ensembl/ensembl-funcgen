-- This one is for general patches
-- Can be applied multiple times without any problems
-- Put dependent patches in later patches

--Some v49 patches which were overwritten 
alter table data_set modify name varchar(100) default NULL;
alter table feature_set modify name varchar(100) default NULL;

-- Alter feature_type and add PSEUDO feature types for regulatory string
alter table feature_type modify `class` enum('Insulator','DNA','Regulatory Feature','Histone','RNA','Polymerase','Transcription Factor','Transcription Factor Complex','Overlap','Regulatory Motif','Region','Enhancer','Expression', 'Pseudo') default NULL;



--New patches - only done on human

--add description to feature_set
alter table feature_set add `description` varchar(80) default NULL,

-- add result_feature column to result_set
--alter table result_set add `result_feature_set` tinyint(1) unsigned default 0;  
--Now done with status

-- add result_feature table
CREATE TABLE `result_feature` (
  `result_feature_id` int(10) unsigned NOT NULL auto_increment,
  `result_set_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) NOT NULL,
  `seq_region_end` int(10) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `window_size` smallint(5) unsigned NOT NULL,
  `score` blob,
  PRIMARY KEY  (`result_feature_id`),
  KEY `set_window_seq_region_idx` (`result_set_id`, `window_size`,`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=50;


--alter table experiment change `date` `created` datetime DEFAULT CURRENT_TIMESTAMP;
--alter table experiment change `date` `created` date DEFAULT CURRENT_DATE;
-- These are invalid defaults?!


--To do
--Add analysis_id to experimental_set, mirror ResultSet unique ley to enable non-unique names with different analysis/cell/feature_types?
--add description to feature_set to enable version annotation for external sets?
--xref_alignment
--This would enable storage of to probe/transcript alignments
-- can we just use identity xref?

--do we need to add array_chip_id to probeset?
--This would only be useful if we want to retrieve all the probesets for an array chip, but not necessarily all of the probes?
--The only way to get all array names/ids, is to retrieve all the probes!!!  Is this any worse than the previous implementation?
--The only answer would be to make probe_set have multiple records for each array_chip_id
--
