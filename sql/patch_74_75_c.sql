/**
@header patch_74_75_c.sql - input_subset.analysis_id_experiment_idx 
@desc   An analysis_id has been added to the input_subset table.
        An experiment_id index has been added.
*/


-- This will mirror the input_set.analysis_id. Consequently, InputSubset has been 
-- changed to inherit from Set, and the feature_type validation in of Set subclass 
-- constructors has been moved to the Set constructor.
-- This work is a prerequisite to the retiring the InputSet class/table.

ALTER TABLE input_subset ADD `analysis_id` smallint(5) unsigned NOT NULL;
ALTER TABLE input_subset ADD KEY analysis_idx (analysis_id);
ALTER TABLE input_subset ADD KEY experiment_idx (experiment_id);


UPDATE input_subset iss, input_set_input_subset isiss, input_set inp set iss.analysis_id=inp.analysis_id where inp.input_set_id=isiss.input_set_id and isiss.input_subset_id=iss.input_subset_id;

analyze table input_subset;
optimize table input_subset;

  
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_c.sql|input_subset.analysis_id_experiment_idx');


