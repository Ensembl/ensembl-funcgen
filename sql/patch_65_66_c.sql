/** 
@header patch_65_66_c.sql - array_chip.name_design
@desc   Increase length of design_id and name fields in array_chip table
*/

ALTER TABLE array_chip MODIFY `design_id` varchar(100) default NULL; 
ALTER TABLE array_chip MODIFY  `name` varchar(100) default NULL;

analyze table array_chip;
optimize table array_chip;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_65_66_c.sql|array_chip.name_design');


