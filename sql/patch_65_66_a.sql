/** 
@header patch_65_66_a.sql - schema version
@desc   Update schema_version in meta table to 66
*/

UPDATE meta SET meta_value='66' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_65_66_a.sql|schema_version');


