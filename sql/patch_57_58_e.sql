# patch_57_58_e.sql
#
# title: probe.name_key
#
# description:
# Add missing name key

ALTER table probe ADD KEY `name_idx` (`name`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_57_58_e.sql|probe.name_key');

