/**
@header patch_61_62_j - meta_key.extend

@desc   Extend the width of the meta_key field in the met table
*/

ALTER TABLE meta MODIFY `meta_key` varchar(46) NOT NULL;


# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_61_62_j.sql|meta_key.extend');
