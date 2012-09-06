/** 
@header patch_68_69_d.sql - xref.id_index_fix
@desc   Update the id_index and info_text definitions to match core schema
*/



ALTER table xref MODIFY `info_text` varchar(255) NOT NULL DEFAULT '';
ALTER table xref drop key id_index;
ALTER table xref ADD UNIQUE KEY `id_index` (`dbprimary_acc`,`external_db_id`,`info_type`,`info_text`,`version`);


OPTIMIZE TABLE xref;
ANALYZE TABLE xref;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_68_69_d.sql|xref.id_index_fix');


