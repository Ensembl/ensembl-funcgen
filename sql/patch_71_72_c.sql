/**
@header patch_71_72_c.sql - supporting_set PK
@desc   Add 'type' to PK of supporting_set PK
*/

ALTER TABLE supporting_set DROP PRIMARY KEY, ADD PRIMARY KEY(`data_set_id`,`supporting_set_id`,`type`);

OPTIMIZE TABLE supporting_set;
ANALYZE TABLE supporting_set;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
 VALUES (NULL, 'patch', 'patch_71_72_c.sql|added_type_to_supporting_set_PK');
