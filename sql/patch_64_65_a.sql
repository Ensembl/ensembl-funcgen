/** 
@header patch_64_65_a.sql - schema version
@desc   Update schema_version in meta table to 65
*/

UPDATE meta SET meta_value='65' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_a.sql|schema_version');


