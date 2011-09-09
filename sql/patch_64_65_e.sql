/** 
@header patch_64_65_e.sql - regulatory_attribute.table_enum
@desc   Add hermaphrodite to the enum for gender
*/

#Need to alter enum here too
ALTER table regulatory_attribute MODIFY `attribute_feature_table` enum('annotated', 'motif') default NULL;

analyze table regulatory_attribute;
optimize table regulatory_attribute;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_e.sql|regulatory_attribute.table_enum');


