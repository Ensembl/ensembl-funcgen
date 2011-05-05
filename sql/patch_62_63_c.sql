/** 
@header patch_62_63_c.sql - binding_matrix.threshold
@desc   Add threshold field to binding_matrix table
*/

ALTER table binding_matrix ADD `threshold` double default NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_62_63_c.sql|binding_matrix.threshold');
