# patch_56_57_f.sql
#
# title: result_feature.scores
#
# description:
# Change result_feature.score to long blob of scores

#Cannot do this as we may have replicate reads
#ALTER table result_feature DROP KEY `set_window_seq_region_idx`;
#ALTER table result_feature ADD UNIQUE KEY `set_window_seq_region_idx` (`result_set_id`, `window_size`,`seq_region_id`,`seq_region_start`);


ALTER TABLE result_feature CHANGE score scores longblob NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_f.sql|result_feature.scores');


 
