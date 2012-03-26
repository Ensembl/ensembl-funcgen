/** 
@header patch_66_67_c.sql - result_feature.remove_partitions
@desc   Remove window_size based partitions as this table should now 
        only hold 0bp data i.e. array data
*/

-- Currently don't use this table but maintain array support for now.

-- Need to conditionally exit if we see window_size !=0
-- Can't do this without defining stored construct (e.g. procedure or function)
-- exit via this construct would leave it in place and requiring tidy up
-- Hence exit, and require manual patch?
-- Rerunning patch after data clean up will remove stored contruct!

--   IF n !=0 THEN select "You have ";
--    END IF;

-- Can't use exit from within procedure or function
-- Will have to wrap up entire patch in function?

DROP PROCEDURE IF EXISTS ApplyPatch;

DELIMITER //

CREATE PROCEDURE ApplyPatch(IN n INT)

  BEGIN
--	DECLARE EXIT HANDLER FOR SQLSTATE '42000'
-- 42000 = procedure does not exist used for non-existant raise_error later
-- Removed this as we actually need to raise and error and exit

    IF n !=0 THEN 
		SELECT "Detecting result_feature records with window_size > 0. Please remove before retrying this patch";
	ELSE
	
		DROP table if exists tmp_result_feature;


		CREATE TABLE `tmp_result_feature` (
			`result_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
			`result_set_id` int(10) unsigned NOT NULL,
			`seq_region_id` int(10) unsigned NOT NULL,
			`seq_region_start` int(10) NOT NULL,
			`seq_region_end` int(10) NOT NULL,
			`seq_region_strand` tinyint(4) NOT NULL,
			`window_size` smallint(5) unsigned NOT NULL,
			`scores` longblob NOT NULL,
		  KEY `result_feature_idx` (`result_feature_id`),
		  KEY `set_window_seq_region_idx` (`result_set_id`,`window_size`,`seq_region_id`,`seq_region_start`)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1;

		INSERT into tmp_result_feature SELECT * from result_feature;
		DROP table result_feature;

	    CREATE table result_feature like tmp_result_feature;
		INSERT into result_feature SELECT * from tmp_result_feature;
		DROP table tmp_result_feature;


		# patch identifier
		INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_66_67_c.sql|result_feature.remove_partitions');


    END IF;

--   IF n !=0 THEN CALL non_existant_procedure;
--   END IF;

  END //

DELIMITER ;

set @num_windows=(select count(distinct window_size) from result_feature where window_size!=0);

call ApplyPatch(@num_windows);

DROP PROCEDURE IF EXISTS ApplyPatch;
