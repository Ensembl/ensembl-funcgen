/**
@header patch_73_74_d.sql - status_name_length
@desc   Increase the length of status_name.name
*/

ALTER TABLE `status_name` MODIFY `name` varchar(60);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_d.sql|status_name_length');


