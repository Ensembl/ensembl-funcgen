/** 
@header patch_63_64_d.sql - experimental_group.url 
@desc   Add URL field to experimental_group
*/

ALTER table experimental_group ADD url varchar(255) default NULL;
ALTER table experimental_group ADD is_project boolean default False;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_63_64_d.sql|experimental_group.url');


