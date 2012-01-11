/** 
@header patch_65_66_d.sql - Reorder an index in unmapped_object.
@desc   The unique_unmapped_obj_idx index in the unmapped_object table is
 		ineffective when querying the table.  Its first part, for example,
 		partly overlaps with the id_idx index.
*/

ALTER TABLE unmapped_object
  DROP INDEX unique_unmapped_obj_idx,
  ADD UNIQUE INDEX unique_unmapped_obj_idx
    (ensembl_id, ensembl_object_type, identifier, unmapped_reason_id,
    parent, external_db_id);

ANALYZE table unmapped_object;
OPTIMIZE table unmapped_object;

# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_65_66_d.sql|unmapped_object.reorder_unmapped_obj_index');
