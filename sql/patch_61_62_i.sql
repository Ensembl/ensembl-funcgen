/**
@header patch_61_62_i - drop_external_db.display_label_linkable

@desc   Remove field display_label_linkable from external_db table as it is obsolete
*/

ALTER TABLE external_db DROP display_label_linkable;

# Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_61_62_i.sql|drop_external_db.display_label_linkable');
