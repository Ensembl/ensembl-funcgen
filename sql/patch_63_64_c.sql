/** 
@header patch_63_64_c.sql - experiment.accession 
@desc   Add archive_id and data_url fields to experiment
*/

ALTER table experiment ADD archive_id varchar(20) default NULL;
ALTER table experiment ADD data_url varchar(255) default NULL;
ALTER table experiment ADD UNIQUE KEY `archive_idx`(`archive_id`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_c.sql|experiment.archive_id_data_url');


