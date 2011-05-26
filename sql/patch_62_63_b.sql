/** 
@header patch_61_62_b.sql - binding_matrix.analysis_id_small_int
@desc   Modify binding_matrix.analysis to smallint to match analysis.analysis_id
*/

ALTER table binding_matrix MODIFY analysis_id smallint(5) unsigned NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_62_63_b.sql|binding_matrix.analysis_id_small_int');


