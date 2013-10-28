/**
@header patch_73_74_c.sql - result_set.replicate
@desc   Add replicate field to result_set
*/
ALTER TABLE `result_set` ADD `replicate` tinyint(3) unsigned NOT NULL;
ALTER TABLE `result_set_input` MODIFY `table_name` enum('experimental_chip','channel','input_set', 'input_subset') DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_c.sql|result_set.replicate');

