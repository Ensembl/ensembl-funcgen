# patch_54_55_d.sql
#
# title: edb.db_name variation fix
#
# description:
# Correct edb db_name for variation db

update external_db set db_name=replace(db_name, '_core_', '_variation_') where db_name like "%_core_Variation";



INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_d.sql|edb.db_name variation fix');

