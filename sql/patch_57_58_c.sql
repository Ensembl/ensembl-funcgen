# patch_57_58_c.sql
#
# title: meta.meta_value_length
#
# description:
# Increase the length of the meta_value field to accomodate long reg strings


ALTER table meta modify `meta_value` varchar(950) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_57_58_c.sql|meta.meta_value_length');


