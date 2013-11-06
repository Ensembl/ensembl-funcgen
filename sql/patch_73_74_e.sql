/**
@header patch_73_74_e.sql - drop_probe_design
desc   Drop the probe_design table

*/

DROP TABLE probe_design;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_73_74_e.sql|drop_probe_design');

