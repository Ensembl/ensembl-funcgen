/** 
@header patch_66_67_d.sql - regulatory_feature.binary_string_500
@desc   Extend binary_string field until we move to binary/blob field
*/



ALTER table regulatory_feature MODIFY binary_string varchar(500) DEFAULT NULL;
analyze table regulatory_feature;
optimize table regulatory_feature;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_66_67_d.sql|regulatory_feature.binary_string_500');


