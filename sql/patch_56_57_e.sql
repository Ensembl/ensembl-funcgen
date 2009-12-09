# patch_56_57_e.sql
#
# title: status.table_name length increase
#
# description:
# Increase max length of status.table_name to accomodate custom data tables


ALTER TABLE status MODIFY `table_name` varchar(32) NOT NULL DEFAULT '';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_e.sql|s.table_name_length');


