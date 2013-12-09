/**
@header patch_74_75_a.sql - schema version
@desc   Update schema_version in meta table to 75
*/

UPDATE meta SET meta_value='75' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_a.sql|schema_version');

