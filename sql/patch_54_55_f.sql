# patch_54_55_f.sql
#
# title: probe_set.name index
#
# description:
# Add missing index

create index name on probe_set(`name`);

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_f.sql|probe_set.name_index');

