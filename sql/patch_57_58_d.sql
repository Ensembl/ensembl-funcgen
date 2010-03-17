# patch_57_58_d.sql
#
# title: meta_key.production_name
#
# description:
# Change meta key species.system_name/ensembl_latin_name to production_name

update meta set meta_key='species.production_name' where meta_key='species.system_name';
update meta set meta_key='species.production_name' where meta_key='species.ensembl_latin_name';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_57_58_d.sql|meta_key.production_name');


