/** 
@header patch_64_65_f.sql - input_set.segmentation_type
@desc   Add segmentation to input_set.type enum
*/

#Need to alter enum here too
ALTER table input_set MODIFY `type` enum('annotated', 'result', 'segmentation') default NULL;

analyze table input_set;
optimize table input_set;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_f.sql|input_set.segmentation_type');


