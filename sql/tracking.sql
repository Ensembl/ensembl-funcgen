-- Reasons for use of separate tracking/stats tables was two fold
-- 1 For visibility (although some of this is now integrated into status)
-- 2 As the data model was originally slightly broken, we needed
--   a way to model many to one relationships between merged input_subsets and 
--   individual records in input_subset_tracking

-- TODO
-- 1 Add support for tracking status history



DROP TABLE IF EXISTS `experiment_tracking`;

CREATE TABLE `experiment_tracking` (
  `experiment_id` int(10) unsigned NOT NULL,
  `notes` TEXT default NULL,
  PRIMARY KEY  (`experiment_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--Slight overkill for one notes field, but keep tables separate for visibility



DROP TABLE IF EXISTS `input_set_tracking`;
/* CREATE TABLE `input_set_tracking` (
  `input_set_id` int(10) unsigned NOT NULL,
  `release_version` varchar(45),
  `is_current` int(1),
  `status` enum('OK','REMOVE','REBUILD') default 'OK',
  PRIMARY KEY  (`input_set_id`),
  UNIQUE KEY `name_cell_type_feature_type_idx` (`experiment_name`,`species`,`cell_type`,`feature_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
*/


DROP TABLE IF EXISTS `experiment_tracking`;

CREATE TABLE `experiment_tracking` (
  `experiment_id` int(10) unsigned NOT NULL,
  `notes` TEXT default NULL,
  PRIMARY KEY  (`experiment_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `result_set_tracking`;

CREATE TABLE `result_set_tracking` (
  `result_set_id`        int(10) unsigned NOT NULL,
  `idr_max_peaks`        mediumint(8) unsigned default NULL,
  `idr_peak_analysis_id` smallint(5) unsigned default NULL,
  PRIMARY KEY  (`result_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- idr stuff in here until we sort the IDRPeaks > DefineOutputSets config linking
-- Also add some summary alignment stats and QC
-- Also separate QC table with full txt file import (for now) and add methods to parse as required
-- this will be cs specific, so probably better to have result/feature_set specific qc tables


-- also add in feature_set_stats here!


DROP TABLE IF EXISTS `input_subset_tracking`;

CREATE TABLE `input_subset_tracking` (
  `input_subset_id` int(10) unsigned NOT NULL,
  `availability_date` date default NULL,
  `download_url` text DEFAULT NULL,
  `download_date` date default NULL,
  `md5sum` varchar(45) default NULL,
  `local_url` text DEFAULT NULL,
  `notes` TEXT default NULL,
  PRIMARY KEY  (`input_subset_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- add notes here as we have lost this from the input_set_tracking table?

-- Change this to current repository...
-- INSERT INTO meta (meta_key, meta_value) VALUES ('fastq_repository','/Your/species/specific/fastq/dir/');
-- INSERT INTO meta (meta_key, meta_value) VALUES ('current_coord_system','GRCh37');

-- Change this to something more appropriate tracking_only? (is_dev is ambiguous)
ALTER TABLE status_name ADD tracking_only TINYINT(1) NOT NULL default 0;

INSERT INTO status_name (name, tracking_only) values ('DOWNLOADED',              1);
INSERT INTO status_name (name, tracking_only) values ('IS_CURRENT',              1);
INSERT INTO status_name (name, tracking_only) values ('ADD_TO_REGULATORY_BUILD', 1);
INSERT INTO status_name (name, tracking_only) values ('IN_REGULATORY_BUILD',     1);
INSERT INTO status_name (name, tracking_only) values ('RELEASED',                1);
INSERT INTO status_name (name, tracking_only) values ('TO_BE_REVOKED',           1);
INSERT INTO status_name (name, tracking_only) values ('REVOKED',                 1);
INSERT INTO status_name (name, tracking_only) values ('TO_BE_REBUILD',           1);
INSERT INTO status_name (name, tracking_only) values ('REBUILT',                 1);
INSERT INTO status_name (name, tracking_only) values ('IN_RELEASE',              1);
INSERT INTO status_name (name, tracking_only) values ('ALIGNING_CONTROL',        1);
INSERT INTO status_name (name, tracking_only) values ('ALIGNED',                 1);
INSERT INTO status_name (name, tracking_only) values ('ALIGNED_CONTROL',                 1);
-- Are these last two ALIGN* states really tracking_only?



-- what about status history?
-- We need to be able to track when a set was released, and whether it was pulled and then re-release
-- Separate table?


--ALTER TABLE `data_set_tracking`     DROP `is_current`;
DROP TABLE IF EXISTS `data_set_tracking`;
CREATE TABLE `data_set_tracking` (
  `data_set_id` int(10) unsigned NOT NULL,
  `release_version` varchar(45) DEFAULT NULL,
  `is_current` int(1) DEFAULT NULL,
  PRIMARY KEY (`data_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- Do we really need data_set_tracking
-- This probably more accruately feature_set_tracking?
-- although we may also want to release stand alone result_sets


--ALTER TABLE `input_set_tracking`    DROP `status`;


--It seems like arrays and optional params or not well supported in functions
--So have 2 functions:
--  1 SumariseRegBuild uses meta string to build summary
--  2 SummariseCellType takes a varchar

-- Could also have SummariseExperiments which would take a FeatureType or a CellType
-- and summarise all the experiments for the given type.


--This is a procedure not a function as we don't return a value
--we simply perform a select

-- We can probably split this up, so we have a *private* procedure to populate
-- cell_type/feature_type records:
--  _UpdateFeatureTypeSummary
--  _UpdateCellTypeSummary
-- As these are private, we can force the use of passing NULL, to mean
-- use the meta entries


-- Or should these have function wrappers? So we don't have to use call?

-- we need separate procedures for Seg and Reg build
-- RegBuild, should have no limitations on feature types
-- Or can we add separate columns to support seg and non-seg ftypes?
-- This might be useful to find celltypes which have a wealth of histones/TFs
-- but a relative paucity of of seg ftypes

DELIMITER //

DROP PROCEDURE IF EXISTS SummariseSegBuild
CREATE PROCEDURE SummariseSegBuild()
BEGIN
/*DECLARE VALUES YOU MAY NEED, EXAMPLE:
  DECLARE NOM_VAR1 DATATYPE [DEFAULT] VALUE;
  */
  --DECLARE CT_IDS VARCHAR DEFAULT NULL;	
  DECLARE @FT_IDS VARCHAR DEFAULT NULL;
  DECLARE @CT_IDS VARCHAR DEFAULT NULL;
  SELECT string INTO @CT_IDS from regbuild_string where name='regbuild.cell_type_ids';
  SELECT string INTO @FT_IDS FROM regbuild_string WHERE name='segbuild.feature_type_ids';		
  
  call _CreateSummaryTable(@FT_IDS, @CT_IDS);
  
END 

-- do we need to allow FT_IDS to be NULL here, in which case we all? For RegBuildSummary?
-- But we really want this to include mandatory segbuild ftypes too even if they are not there.
-- This is just a union, no?


CREATE PROCEDURE _CreateSummaryTable(IN FT_IDS VARCHAR(1000), IN CT_IDS VARCHAR(1000))
BEGIN
	--IF CT_IDS is NULL THEN
    	 
	--END IF;
	
	--FT_IDS Will never be NULL As this is always known by the caller
	
	DROP table if exists `progress_summary`;
	CREATE TABLE `progress_summary` (
	 `cell_type_id`          int(10) unsigned NOT NULL,
	 `feature_type_id`       int(10) unsigned NOT NULL,
	 `experiments`           int(10) DEFAULT NULL,
	 `aligned_result_sets`   int(10) DEFAULT NULL,
	 `imported_feature_sets` int(10) DEFAULT NULL,
	 PRIMARY KEY (`cell_type_id`),
	 UNIQUE KEY `cell_feature_type_idx` (`cell_type_id`, `feature_type_id`)
  ) ENGINE=MyISAM;
	
	
	--This will create enpty records for each ctype/ftype combination to prevent having
	--to do problematic nested left/right joins.
	--The IN will cast the scalar value of CT_IDS as an int and crop it at the first comma
	--If MySQL supported arrays, then this would be possible, and would even use and index
	--but we are limited to FIND_IN_SET here, which is just a string operation to return
	--the index in the string of comma separated values.
	
	INSERT INTO progress_summary(cell_type_id, feature_type_id) 
		SELECT ct.cell_type_id, ft,feature_type_id from cell_type ct, feature_type ft 
			WHERE FIND_IN_SET(ct.cell_type_id, CT_IDS) and FIND_IN_SET(ft.feature_type_id, FT_IDS);
			
	--All Experiments
	UPDATE progress_summary ps, (select cell_type_id, feature_type_id, count(*) as cnt from experiment 
								 group by cell_type_id, feature_type_id) e 
	 SET ps.experiments=e.cnt WHERE FIND_IN_SET(e.cell_type_id, CT_IDS) AND FIND_IN_SET(e.feature_type_id, FT_IDS);
	
	--All ResultSets 
	--do we want to split this into all and aligned?
	--We have a join to experiment here to make handle the IDR replicates i.e. we only wnat to count the distinct experiment
	--we could do this by counting the distinct input_subset.experiment_id where is_control=0
	UPDATE progress_summary ps, 
	 (select cell_type_id, feature_type_id, (IFNULL( COUNT(exp_id) , 0 )) cnt from 
	   (select distinct(iss.experiment_id) as exp_id, rs.feature_type_id, rs.cell_type_id 
	    from input_subset iss join result_set_input rsi on iss.input_subset_id=rsi.table_id 
	    join result_set rs using(result_set_id) where rsi.table_name='input_subset') 
    rse group by cell_type_id, feature_type_id)
	SET ps.experiments=rs.cnt WHERE FIND_IN_SET(rs.cell_type_id, CT_IDS) AND FIND_IN_SET(rs.feature_type_id, FT_IDS);
	 
	 --select cell_type_id, feature_type_id, (IFNULL( COUNT(exp_name) , 0 )) cnt from (select distinct(e.name) as exp_name, e.feature_type_id, rs.cell_type_id from experiment e left join result_set rs on rs.name like concat(e.name, '%') where rs.name is not NULL) rse group by cell_type_id, feature_type_id;
	 
	--How do we then dynamically create the output table query?
	--This is easy as we want to grep by ct ft, and then simply have the counts as the field headers
	
	SELECT ct.name, ft.name, ps.experiments, ps.aligned_result_sets, ps.imported_feature_sets 
	 FROM feature_type ft 
	 JOIN progress_summary ps USING(feature_fype_id)
	 JOIN cell_type ct    using (cell_type_id);
	
	
END 



//
DELIMITER;

