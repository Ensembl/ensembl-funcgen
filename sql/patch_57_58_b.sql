# patch_57_58_b.sql
#
# title: regulatory_feature.cell_type
#
# description:
# Add cell_type_id and associated indexes to regulatory_feature


ALTER table regulatory_feature ADD column `cell_type_id` int(10) unsigned DEFAULT NULL AFTER feature_type_id;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_57_58_b.sql|regulatory_feature.cell_type');


