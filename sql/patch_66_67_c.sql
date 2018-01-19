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
		
		-- removing window_size, but leaving scores as blob as
		-- conversion of data connot be done in this patch
		-- and this table will most likely be removed shortly

		CREATE TABLE `tmp_result_feature` (
			`result_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
			`result_set_id` int(10) unsigned NOT NULL,
			`seq_region_id` int(10) unsigned NOT NULL,
			`seq_region_start` int(10) NOT NULL,
			`seq_region_end` int(10) NOT NULL,
			`seq_region_strand` tinyint(4) NOT NULL,
			`scores` longblob NOT NULL,
		  PRIMARY KEY `result_feature_idx` (`result_feature_id`),
		  KEY `set_seq_region_idx` (`result_set_id`,`seq_region_id`,`seq_region_start`)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1;

		INSERT into tmp_result_feature 
			SELECT result_feature_id, result_set_id, seq_region_id, seq_region_start,
					 seq_region_end, seq_region_strand, scores
			FROM result_feature;
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
