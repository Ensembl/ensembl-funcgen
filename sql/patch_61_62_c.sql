/** 
@header patch_61_62_c.sql - feature_type.sequence_ontology
@desc   Add so_accession and so_name and fields
*/


#Or generic associated_ontology table?


ALTER table feature_type ADD `so_accession` varchar(64) DEFAULT NULL;
ALTER table feature_type ADD `so_name` varchar(255) DEFAULT NULL;

ALTER table feature_type ADD KEY `so_accession_idx` (`so_accession`);

analyze table feature_type;
optimize table feature_type;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_61_62_c.sql|feature_type.sequence_ontology');
