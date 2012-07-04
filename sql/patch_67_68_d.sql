/** 
@header patch_67_68_d.sql - feature_set.input_set_id
@desc   Add input_set_id in place of experiment_id
*/


ALTER TABLE feature_set ADD `input_set_id` int(10) unsigned DEFAULT NULL;


UPDATE feature_set fs, input_set inp set fs.input_set_id=inp.input_set_id where fs.name like concat(inp.name, '%');

ALTER TABLE feature_set DROP `experiment_id`;
-- This drops the experiment_id index

-- Not replacing with input_set_id index as we currently don't query with an InputSet

OPTIMIZE TABLE feature_set;
ANALYZE TABLE feature_set;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_67_68_d.sql|feature_set.input_set_id');


