# patch_56_57_k.sql
#
# title: meta_species.system_name
#
# description:
# Change species.ensembl_latin_name to system_name


update meta set meta_key='species.system_name' where meta_key='species.ensembl_latin_name';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_k.sql|meta_species.system_name');


 
