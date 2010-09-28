# patch_59_60_h.sql
#
# Title: af_amf.index_tweaks
#   
#
# Description:
#   Alter annotated_feature and motif_feature indexes to increase performance


ALTER TABLE associated_motif_feature ADD KEY `motif_feature_idx` (`motif_feature_id`);
OPTIMIZE TABLE associated_motif_feature;

ALTER TABLE annotated_feature ADD UNIQUE KEY `seq_region_feature_set_idx` (`seq_region_id`,`seq_region_start`,`feature_set_id`);
ALTER TABLE annotated_feature DROP KEY seq_region_idx; 
OPTIMIZE TABLE annotated_feature;

# Insert patch identifier.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_59_60_h.sql|af_amf.index_tweaks');
