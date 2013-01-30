/** 
@header patch_70_71_c.sql - design_table_removal
@desc   Remove unused design tables
*/

DROP TABLE design_type;
DROP TABLE experimental_design;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_70_71_c.sql|design_table_removal');


