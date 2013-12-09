/**
@header patch_74_75_b.sql - result_set.name_unique
@desc   Patch to make name field of result_set table unique, by appending analysis
        logic_name.
*/

UPDATE result_set rs, analysis a set rs.name=concat(rs.name, '_', a.logic_name) where rs.analysis_id=a.analysis_id;

ALTER TABLE result_set DROP INDEX unique_dx;
ALTER TABLE result_set ADD UNIQUE KEY name_idx (name);

-- Now add back in other keys.

ALTER TABLE result_set ADD KEY cell_type_idx (cell_type_id);
ALTER TABLE result_set ADD KEY feature_type_idx (feature_type_id);
ALTER TABLE result_set ADD KEY analysis_idx (analysis_id);
ALTER TABLE result_set ADD KEY feature_class_idx (feature_class);


analyze table result_set;
optimize table result_set;

  
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_b.sql|result_set.name_unique');


