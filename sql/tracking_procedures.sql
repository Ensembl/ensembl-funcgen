

-- It seems like arrays and optional params or not well supported in functions
-- So have 2 functions:
--  1 SumariseRegBuild uses meta string to build summary
--  2 SummariseCellType takes a varchar

-- Could also have SummariseExperiments which would take a FeatureType or a CellType
-- and summarise all the experiments for the given type.


-- This is a procedure not a function as we don't return a value
-- we simply perform a select

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

-- to do, turn this in to a func, so we don't have to use 'call'

DROP PROCEDURE IF EXISTS SummariseSegBuild;
DELIMITER //
CREATE PROCEDURE SummariseSegBuild()
SummariseSegBuild: BEGIN
  DECLARE FT_IDS, CT_IDS VARCHAR(1000) DEFAULT NULL;
  SET CT_IDS=GetRegBuildCellTypeIDs();
  SET FT_IDS=GetSegBuildFeatureTypeIDs();
  
  
  IF CT_IDS is NULL THEN
    SELECT 'The regbuild_string regbuild.cell_type_ids is not defined' as '';    
    LEAVE SummariseSegBuild;
  END IF;
  
  IF FT_IDS is NULL THEN
    SELECT 'The regbuild_string segbuild.feature_type_ids is not defined' as '';
    LEAVE SummariseSegBuild;
  END IF;
  
  
  call _CreateSummaryTable(FT_IDS, CT_IDS);  
END 
//
DELIMITER ;

-- do we need to allow FT_IDS to be NULL here, in which case we all? For RegBuildSummary?
-- But we really want this to include mandatory segbuild ftypes too even if they are not there.
-- This is just a union, no?


DROP FUNCTION IF EXISTS GetSegBuildFeatureTypeIDs;
DELIMITER //
CREATE FUNCTION GetSegBuildFeatureTypeIDs() returns VARCHAR(1000)
BEGIN
	-- apparently can't use select FT_IDS := string
	-- unless we are using a @FT_IDS style variable (which is a user session var)
	DECLARE FT_IDS VARCHAR(1000) DEFAULT NULL;
	SELECT string into FT_IDS FROM regbuild_string WHERE name='segmentation.feature_type_ids';
	RETURN FT_IDS;
END 
//
DELIMITER ;

DROP FUNCTION IF EXISTS GetRegBuildCellTypeIDs;
DELIMITER //
CREATE FUNCTION GetRegBuildCellTypeIDs() returns VARCHAR(1000)
BEGIN
  DECLARE CT_IDS VARCHAR(1000) DEFAULT NULL;
  SELECT string into CT_IDS FROM regbuild_string WHERE name='regbuild.cell_type_ids';
  RETURN CT_IDS;
END 
//
DELIMITER ;


DROP table if exists progress_summary;
  CREATE TABLE `progress_summary` (
   `cell_type_id`          int(10) unsigned NOT NULL,
   `feature_type_id`       int(10) unsigned NOT NULL,
   `experiments`           int(10) DEFAULT NULL,
   `imported_feature_sets` int(10) DEFAULT NULL,
   `imported_result_sets`  int(10) DEFAULT NULL,
   `aligned_result_sets`   int(10) DEFAULT NULL,
   PRIMARY KEY (`cell_type_id`, `feature_type_id`)
  ) ENGINE=MyISAM;

CREATE OR REPLACE VIEW progress_view AS
   SELECT ct.name as CellType, ft.name as FeatureType, ps.experiments, (ps.experiments - ps.imported_result_sets) as unprocessed_experiments, ps.imported_feature_sets, ps.aligned_result_sets, ps.imported_result_sets
   FROM feature_type ft 
   JOIN progress_summary ps USING(feature_type_id)
   JOIN cell_type ct    using (cell_type_id);

/* This requires a procedure using separate updates as the 
 * left join between multiple tables broke the query
 * and MySQL does not support right joins
 * i.e. we were losing the NULL rows
 */

DROP PROCEDURE IF EXISTS _CreateSummaryTable;
DELIMITER //
CREATE PROCEDURE _CreateSummaryTable(IN FT_IDS VARCHAR(1000), IN CT_IDS VARCHAR(1000))
_CreateSummaryTable: BEGIN
	-- IF CT_IDS is NULL THEN
	-- END IF;
	-- FT_IDS Will never be NULL As this is always known by the caller
	
	-- SELECT '_CreateSummaryTable FT_IDS', FT_IDS;
	-- SELECT '_CreateSummaryTable CT_IDS', CT_IDS;

	-- This will create enpty records for each ctype/ftype combination to prevent having
	-- to do problematic nested left/right joins.
	-- The IN will cast the scalar value of CT_IDS as an int and crop it at the first comma
	-- If MySQL supported arrays, then this would be possible, and would even use and index
	-- but we are limited to FIND_IN_SET here, which is just a string operation to return
	-- the index in the string of comma separated values.
	
	truncate progress_summary;
	INSERT INTO progress_summary(cell_type_id, feature_type_id) 
		SELECT ct.cell_type_id, ft.feature_type_id from cell_type ct, feature_type ft 
			WHERE FIND_IN_SET(ct.cell_type_id, CT_IDS) and FIND_IN_SET(ft.feature_type_id, FT_IDS) order by ct.name, ft.name;
			
			
	/*	SELECT ps.cell_type_id, ps.feature_type_id, e.cnt from progress_summary ps, (select cell_type_id, feature_type_id, count(*) as cnt from experiment 
          WHERE FIND_IN_SET(cell_type_id, CT_IDS) AND FIND_IN_SET(feature_type_id, FT_IDS) group by cell_type_id, feature_type_id) e where e.cell_type_id=ps.cell_type_id 
          and ps.feature_type_id=e.feature_type_id; */
          
	-- All Experiments
  UPDATE progress_summary ps, 
	       (select cell_type_id, feature_type_id, count(*) as cnt from experiment 
          WHERE FIND_IN_SET(cell_type_id, CT_IDS) AND FIND_IN_SET(feature_type_id, FT_IDS) group by cell_type_id, feature_type_id) e 
  SET ps.experiments=e.cnt WHERE ps.cell_type_id=e.cell_type_id and ps.feature_type_id=e.feature_type_id;
  UPDATE progress_summary set experiments=0 where experiments is NULL;
  
  
  
	-- All ResultSets 
	-- do we want to split this into all and aligned?
	-- We have a join to experiment here to make handle the IDR replicates i.e. we only wnat to count the distinct experiment
	-- we could do this by counting the distinct input_subset.experiment_id where is_control=0
	UPDATE progress_summary ps, 
	 (select rse.cell_type_id, rse.feature_type_id, (IFNULL( COUNT(distinct(rse.exp_id)) , 0 )) cnt from 
	   (select iss.experiment_id as exp_id, rs.feature_type_id, rs.cell_type_id 
	    from input_subset iss join result_set_input rsi on iss.input_subset_id=rsi.table_id 
	    join result_set rs using(result_set_id) where rsi.table_name='input_subset' AND iss.is_control != 1 AND
	    FIND_IN_SET(rs.cell_type_id, CT_IDS) AND FIND_IN_SET(rs.feature_type_id, FT_IDS)) 
    rse group by cell_type_id, feature_type_id) rs1
	SET ps.imported_result_sets=rs1.cnt WHERE ps.cell_type_id=rs1.cell_type_id and ps.feature_type_id=rs1.feature_type_id;
	UPDATE progress_summary set imported_result_sets=0 where imported_result_sets is NULL;
	 

	-- tricky to handle aligned sets here, as we need to handle replicate sets
	

  -- For the feature_set counts, we still need to link back to the experiment table somehow
  -- to ensure that we are not counting any data which is not associated with an experiment
  -- Unlikely, but the schema/API do support it	
	
	
	
	
	-- How do we then dynamically create the output table query?
	-- This is easy as we want to grep by ct ft, and then simply have the counts as the field headers
	
	-- seemingly can do this just 'if not exists'
	-- meaning this will be redefined ever time this prodecude is called
	-- actually can I create it anyway first?
	
	

	 SELECT * from progress_view;  
	 
	
END 
//
DELIMITER ;



-- we need a procedure to re-create this if the CT/FT_IDS ever change
-- or can we use functions in view to make them dynamic?

-- update these to show the relevant controls in a 3rd column?
-- would be nice to get 'with rollup' in here, but that needs a group by.
-- can we add that artificially?
-- This does not work select * from progress_view group by CellType, FeatureType with rollup;

-- we really need a view which shows whether subsets are downloaded or not

/* Show just the data required for the segmentation data
 * Use this to kick off prioritised jobs based on a group_concat query:
 *  select group_concat('"', experiment, '"') from seg_exp_view 
 *  where experiment like "CTYPE_%_EXP_GROUP_" and result_set is NULL;
 * Or use this wildcard string directly in SeedPipeline. This may identify
 * additional non-segmentation experiments, so it's a good idea to check
 * the reg_exp_view, and use the group_concat approach if there are 
 * already some processed experiments. 
 */


/* Append suggested control sets based on a self join using the experimental_group, cell_type and WCE feature type
 * This needs to be subgroup!
 * This will work once we have patched the experimental_group table
 * To represent subgroups as groups, but we should preferentially use the project 
 * which will need changing from a boolean to a varchar 
 */

DELIMITER //
CREATE OR REPLACE VIEW seg_exp_view AS
  SELECT e.experiment_id, e.name as experiment, ft.class, rs.name as result_set, et.notes 
    FROM feature_type ft, experiment e LEFT JOIN result_set rs on rs.name like concat(e.name, '%') 
    LEFT JOIN experiment_tracking et on e.experiment_id=et.experiment_id
    WHERE FIND_IN_SET(e.cell_type_id, GetRegBuildCellTypeIDs()) 
    AND FIND_IN_SET(e.feature_type_id, GetSegBuildFeatureTypeIDs()) 
    AND e.feature_type_id=ft.feature_type_id ORDER by e.name;
//
DELIMITER ;

/* Shows all the data viable for the regulatory build i.e. also including the 
 * histones and TFs not required for the segmentation
 */

/* TODO implement this without hardcoding for WCE, i.e. use is_control somehow
   We really need to have is_control at the experiment or feature_type level too.
*/
/*
DELIMITER //
CREATE OR REPLACE VIEW reg_exp_view AS
  SELECT e.experiment_id, e.name as experiment, ft.class, rs.name as result_set, et.notes, group_concat(e1.name) as controls
  FROM feature_type ft, experiment e LEFT JOIN result_set rs on rs.name like concat(e.name, '%') 
    LEFT JOIN experiment_tracking et on e.experiment_id=et.experiment_id
    LEFT JOIN experiment e1 on (e.experimental_group_id=e1.experimental_group_id
      AND e.cell_type_id=e1.cell_type_id and e1.name like "%\_WCE\_%")
    WHERE FIND_IN_SET(e.cell_type_id, GetRegBuildCellTypeIDs()) 
    AND e.feature_type_id=ft.feature_type_id GROUP by e.experiment_id, rs.name ORDER by e.name;
//
DELIMITER ;
*/

DELIMITER //
CREATE OR REPLACE VIEW reg_exp_view AS
  SELECT e.experiment_id, e.name as experiment, ft.class, rs.name as result_set, et.notes
  FROM feature_type ft, experiment e LEFT JOIN result_set rs on rs.name like concat(e.name, '%') 
    LEFT JOIN experiment_tracking et on e.experiment_id=et.experiment_id
    WHERE FIND_IN_SET(e.cell_type_id, GetRegBuildCellTypeIDs()) 
    AND e.feature_type_id=ft.feature_type_id ORDER by e.name;
//
DELIMITER ;



DELIMITER //
CREATE OR REPLACE VIEW exp_subset_view AS
  SELECT e.name as experiment, iss.input_subset_id, iss.name, iss.replicate, iss.is_control, isst.availability_date, isst.md5sum, group_concat(sn.name) as states, isst.notes, isst.download_url, isst.local_url  
    FROM experiment e 
    LEFT JOIN input_subset iss using(experiment_id) 
    LEFT JOIN input_subset_tracking isst using(input_subset_id)
    LEFT JOIN status s on iss.input_subset_id=s.table_id
    LEFT JOIN status_name sn on s.status_name_id=sn.status_name_id
    WHERE s.table_name='input_subset'
    GROUP BY iss.input_subset_id;
//
DELIMITER ;


DELIMITER //
CREATE OR REPLACE VIEW result_subset_view AS
  SELECT rs.result_set_id, rs.name AS result_set, iss.input_subset_id, iss.name AS input_subset, 
         iss.replicate, iss.is_control
  FROM   result_set rs 
  JOIN   result_set_input rsi ON rs.result_set_id = rsi.result_set_id
  JOIN   input_subset iss     ON iss.input_subset_id = rsi.table_id
  WHERE  rsi.table_name = 'input_subset';
//
DELIMITER ;

DELIMITER //
CREATE OR REPLACE VIEW result_replicate_subset_view AS
  SELECT rs.result_set_id, rs.name as result_set, iss.name as input_subset, iss.replicate, 
          iss.is_control, concat(replace(rs.name, 'TR', ''), '.bam') as align_file     
  FROM result_set rs      
  JOIN result_set_input rsi USING(result_set_id)      
  JOIN input_subset iss on iss.input_subset_id = rsi.table_id       
  WHERE rsi.table_name='input_subset' and rs.name rlike ".*_TR[0-9]+$";
//
DELIMITER ;

/* we don't really need result_merged_subset_view as the mapping between bams and input fastqs is 
 * only really needed for the segmentation i.e. individual replicates bams
 * this would involved a replcace(group_concat(replicate), ',', '_') and an order on the 
 * replicate aswell as the necessary group keys. Would also need to group concat the input_subset names
 * so not ideal. 
 */

DELIMITER //
CREATE OR REPLACE VIEW exp_status_view AS
  SELECT e.name as experiment, sn.name
    FROM experiment e 
    LEFT JOIN status s on s.table_id=e.experiment_id 
    LEFT JOIN status_name sn using(status_name_id)
    WHERE s.table_name='experiment' ORDER BY e.name;
//
DELIMITER ;

