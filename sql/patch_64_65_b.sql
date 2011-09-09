/** 
@header patch_64_65_b.sql - feature_type.analysis_id
@desc   Add analysis_id field and add to unique key. To support SegmentationFeature states
*/

ALTER table feature_type ADD `analysis_id` smallint(5) unsigned default NULL AFTER class;
ALTER table feature_type DROP KEY `name_class_idx`;
ALTER table feature_type ADD UNIQUE KEY `name_class_analysis_idx` (`name`,`class`, `analysis_id`);


#Need to alter enum here too
ALTER table feature_type MODIFY  `class` enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Regulatory Motif',  'Enhancer', 'Expression', 'Pseudo', 'Open Chromatin', 'Search Region', 'Association Locus', 'Segmentation State') default NULL;

analyze table feature_type;
optimize table feature_type;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_64_65_b.sql|feature_type.analysis_id');


