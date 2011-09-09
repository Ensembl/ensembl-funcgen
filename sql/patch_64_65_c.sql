/** 
@header patch_64_65_c.sql - cell_type.gender_hermaphrodite
@desc   Add hermaphrodite to the enum for gender
*/

#Need to alter enum here too
ALTER table cell_type MODIFY   `gender` enum('male','female', 'hermaphrodite') DEFAULT NULL;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_c.sql|cell_type.gender_hermaphrodite');


