# patch_60_61_b.sql
#
# title:  add binding_matrix.analysis_id and drop binding_matrix.type
#
# description:
# Add analysis_id field to probe table and drop type field


ALTER table binding_matrix ADD `analysis_id` int(10) unsigned NOT NULL;
ALTER table binding_matrix DROP `type`;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_60_61_b.sql|binding_matrix.analysis_id');


