# patch_58_59_e.sql
#
# title: add meta schema_type
#
# description:
# Add meta entry denoting funcgen schema_type

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'schema_type', 'funcgen');
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_58_59_e.sql|meta.schema_type');


