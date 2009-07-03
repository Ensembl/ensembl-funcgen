# patch_55_56_a.sql
#
# title: probe_set.name length
#
# description:
# Extend probe_set.name field length

alter table probe_set modify `name` varchar(100) NOT NULL;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_55_56_a.sql|probe_set.name_length');

