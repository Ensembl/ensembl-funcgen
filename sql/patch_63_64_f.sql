/** 
@header patch_63_64_f.sql - feature_set.experiment_id
@desc   Add experimenta_id to feature_set to enable easy recall of experiment
*/

ALTER table feature_set ADD experiment_id int(10) unsigned default '0';
ALTER table feature_set ADD KEY `experiment_idx` (`experiment_id`);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_f.sql|feature_set.experiment_id');


