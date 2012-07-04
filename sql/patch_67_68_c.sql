/** 
@header patch_67_68_c.sql - input_set_subset.replicate_is_control
@desc   Add replicate field to input_set and input_subset
*/



ALTER TABLE input_subset ADD `replicate` tinyint(3) unsigned NOT NULL;
ALTER TABLE input_subset ADD `is_control` tinyint(3) unsigned NOT NULL;
ALTER TABLE input_set ADD `replicate` tinyint(3) unsigned NOT NULL;

-- If input_subset.replicate=0 and input_subset.is_control=0
-- then current subset file represents a replicate pool
-- with or without a control file

-- If input_subset.replicate=0 and input_subset.is_control=1
-- subset is a pool of controls

-- If input_subset.replicate>0 and input_subset.is_control=1
-- subset is a single control file


OPTIMIZE TABLE input_subset;
ANALYZE TABLE input_set;
OPTIMIZE TABLE input_set;
ANALYZE TABLE input_subset;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_67_68_c.sql|input_set_subset.replicate_is_control');


