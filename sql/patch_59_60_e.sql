# patch_59_60_e.sql
#
# title: regulatory_attribute.allow_motif_table
#
# description:
# Add motif to attribute_feature_table enum

alter table regulatory_attribute MODIFY `attribute_feature_table` enum('annotated','external', 'motif') NOT NULL DEFAULT 'annotated';


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_59_60_e.sql|regulatory_attribute.allow_motif_table');


