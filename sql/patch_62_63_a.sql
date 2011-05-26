/** 
@header patch_62_63_a.sql - schema version
@desc   Update schema_version in meta table to 63
*/

UPDATE meta SET meta_value='63' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_62_63_a.sql|schema_version');


