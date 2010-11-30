# patch_60_61_c.sql
#
# title: probe name extension
#
# description:
# Extend probe name to 100 characters

alter table probe modify `name` varchar(100) NOT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_60_61_c.sql|probe.name_alter');


