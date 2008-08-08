# patch_50_51_d.sql
#
# Title: mirror multispecies changes in core meta
#
# Description:
#   Add a species_id column to the meta, coord_system left for now
#   as we have our own eFG CoordSystemAdaptor

-- Add the new species_id column after meta_id
ALTER TABLE meta ADD COLUMN
 species_id INT UNSIGNED DEFAULT 1 -- Default species_id is 1
                                   -- NULL means "not species specific"
 AFTER meta_id;

-- Redo the indexes on the meta table
ALTER TABLE meta DROP INDEX key_value;
ALTER TABLE meta DROP INDEX meta_key_index;
ALTER TABLE meta DROP INDEX meta_value_index;

ALTER TABLE meta
 ADD UNIQUE INDEX species_key_value_idx (species_id, meta_key, meta_value);
ALTER TABLE meta
 ADD INDEX species_value_idx (species_id, meta_value);

-- Optimize the modified tables
OPTIMIZE TABLE meta;


UPDATE  meta SET species_id = NULL WHERE meta_key IN ('patch', 'schema_version');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_b.sql|multispecies');
