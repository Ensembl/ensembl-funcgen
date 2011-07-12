/** 
@header patch_63_64_b.sql - cell_type.efo_id
@desc   Add Experimental Factor Ontology ID field to cell_type
*/

ALTER table cell_type ADD efo_id varchar(20) default NULL;
ALTER table cell_type ADD UNIQUE KEY `efo_idx`(`efo_id`);

# Is actually only 11 chars long


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_b.sql|cell_type.efo_id');


