/**
@header patch_73_74_c.sql - result_set.replicate
@desc   Add replicate field to result_set
*/
ALTER TABLE `result_set` ADD `replicate` tinyint(3) unsigned NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_c.sql|result_set.replicate');

