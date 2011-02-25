/** 
@header patch_61_62_h - external_db.db_name_idx 

@desc   Add unique index to db_name field in external_db table.
*/

CREATE INDEX db_name_idx ON external_db(db_name);

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_61_62_h.sql|external_db.db_name_idx');
