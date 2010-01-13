# patch_56_57_b.sql
#
# title: unmapped_object object_type key
#
# description:
# Add a ensembl_object_type/id key to the unmapped_object table

#ALTER TABLE unmapped_object DROP KEY `object_type_idx`;

ALTER TABLE unmapped_object ADD KEY `object_type_idx` (`ensembl_id`, `ensembl_object_type`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_b.sql|uo.object_type_id_key');


