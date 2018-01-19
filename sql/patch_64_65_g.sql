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
@header patch_64_65_g.sql - MAX_ROWS_AVG_ROW_LENGTH_clean_up
@desc   Correcting/Updating some MAX_ROWS and AVG_ROW_LENGTH table options
*/


-- remove AVG_ROW_LENGTH spec from some DBs

DROP TABLE IF EXISTS `regulatory_attribute_tmp`;
CREATE TABLE `regulatory_attribute_tmp` (
  `regulatory_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated', 'motif') default NULL,
  PRIMARY KEY  (`regulatory_feature_id`, `attribute_feature_table`, `attribute_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

insert into regulatory_attribute_tmp select * from regulatory_attribute;

drop table  regulatory_attribute;

CREATE TABLE `regulatory_attribute` (
  `regulatory_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated', 'motif') default NULL,
  PRIMARY KEY  (`regulatory_feature_id`, `attribute_feature_table`, `attribute_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

insert into regulatory_attribute select * from regulatory_attribute_tmp;



drop table  regulatory_attribute_tmp;
optimize table regulatory_attribute;
analyze table regulatory_attribute;



-- ALTER table annotated_feature AVG_ROW_LENGTH=80;
-- Or just ignore this and drop display_label and AVG_ROW_LENGTH in v66
-- This field isn't really used at all

ALTER table regulatory_feature MAX_ROWS=100000000;
optimize table regulatory_feature;
analyze table regulatory_feature;


ALTER table object_xref AVG_ROW_LENGTH=40;
ALTER table object_xref MAX_ROWS=100000000;
optimize table object_xref;
analyze table object_xref;


-- actually avg length based on table status info is 41, we'll let it off here

-- Don't have MAX_ROWS with object_xref?
-- So what use it AVG_ROW_LENGTH?
-- Actual entries 21405370
-- Proposed MAX_ROWS 100000000;
-- Would require ~5* increase in oxs. Fine




-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_g.sql|MAX_ROWS_AVG_ROW_LENGTH_clean_up');






-- Do this for v66



-- MAX_ROWS/AVG_ROW_LENGTH for identity_xref and others?
-- Setting to 40 and 100000000
-- Changes the Max_data_length
-- from 281474976710655
-- to   4294967295
-- Use 5* or 100000000 as MAX_ROWS


-- xref
-- For human data suggests AVG_ROW_LENGTH=53;
-- These will keep increasing unless we clear out old xref records with no object_xrefs
-- Does FuncgenForeignKeys HC catch this?
-- ALTER table xref AVG_ROW_LENGTH=53;
-- ALTER table xref MAX_ROWS=100000000;
-- optimize table xref;
-- analyze table xref;


-- unmapped_object
-- 691766060/17107478 = 40.43
-- ALTER table unmapped_object AVG_ROW_LENGTH=40;
-- ALTER table unmapped_object MAX_ROWS=100000000;
-- optimize table unmapped_object;
-- analyze table unmapped_object;

