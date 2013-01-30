/** 
@header patch_70_71_b.sql - analysis_key_clean
@desc   Remove duplicate logic_name key
*/

ALTER TABLE analysis DROP INDEX `logic_name_idx`;
ALTER TABLE analysis DROP INDEX `logic_name`;

ALTER TABLE analysis ADD UNIQUE INDEX logic_name_idx (logic_name);

ANALYZE table analysis;
OPTIMIZE table analysis;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_b.sql|analysis_key_clean');


