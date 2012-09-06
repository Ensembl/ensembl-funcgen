/** 
@header patch_68_69_c.sql - regbuild_string.species_id_not_null
@desc   Change the regbuild_string.species_id to be not null smallint
*/

-- Note species_id is simply an int in meta, but probably needs changing

ALTER table regbuild_string MODIFY species_id smallint(5) unsigned NOT NULL DEFAULT 1;

OPTIMIZE TABLE regbuild_string;
ANALYZE TABLE regbuild_string;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_68_69_c.sql|regbuild_string.species_id_not_null');


