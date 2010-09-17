# patch_59_60_f.sql
#
# title: bm.frequencies_pf_index_mods
#
# description:
# Make some minor modifications to the binding_matric frequences column
# and probe_feature indexes to increase performance

alter table binding_matrix modify `frequencies` varchar(1000) NOT NULL;

alter table probe_feature add index `seq_region_probe_probe_feature_idx` (`seq_region_id`,`seq_region_start`, `seq_region_end`, `probe_id`, `probe_feature_id`); 
alter table probe_feature drop index seq_region_idx; 


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_59_60_f.sql|bm.frequencies_pf.index_mods');


