/** 
@header patch_63_64_e.sql - object_xref.ensembl_object_type 
@desc   Add Experiment to the type of ensembl_object_type field
*/

ALTER table object_xref MODIFY ensembl_object_type ENUM('Experiment', 'RegulatoryFeature', 'ExternalFeature', 'AnnotatedFeature', 'FeatureType', 'ProbeSet', 'Probe', 'ProbeFeature') not NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_e.sql|object_xref.ensembl_object_type');


