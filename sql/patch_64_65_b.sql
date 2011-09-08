/** 
@header patch_64_65_b.sql - feature_type.analysis_id
@desc   Add analysis_id field and add to unique key. To support SegmentationFeature states
*/

ALTER table feature_type ADD `analysis_id` smallint(5) unsigned default NULL AFTER class;
ALTER table feature_type DROP KEY `name_class_idx`;
ALTER table feature_type ADD UNIQUE KEY `name_class_analysis_idx` (`name`,`class`, `analysis_id`);

# Is actually only 11 chars long


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_b.sql|feature_type.analysis_id');


