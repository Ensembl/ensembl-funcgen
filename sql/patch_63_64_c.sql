/** 
@header patch_63_64_c.sql - experiment.accession 
@desc   Add Accession ID and Data URL fields to experiment
*/

ALTER table experiment ADD accession_id varchar(20) default NULL;
ALTER table experiment ADD data_url varchar(255) default NULL;
ALTER table experiment ADD UNIQUE KEY `accession_idx`(`accession_id`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_c.sql|experiment.accession');


