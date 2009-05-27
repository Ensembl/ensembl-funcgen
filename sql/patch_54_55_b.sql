# patch_54_55_a.sql
#
# title: FK consistency tweaks
#
# description:
# Some general schema tweaks to bring var types in line with core and internally consistent
# To allow valid foreign keys when using innodb

-- Change analysis_id to smallint(5)
alter table analysis modify `analysis_id` smallint(5) unsigned NOT NULL auto_increment;
alter table probe_feature modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table feature_set modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table result_set modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table probe_design modify `analysis_id` smallint(5) unsigned NOT NULL;

-- Change s.status_name_id to unsigned
alter table status modify `status_name_id` int(10) unsigned NOT NULL;

--  Change c.coord_system_id to unsigned
alter table coord_system modify `coord_system_id` int(10) unsigned NOT NULL auto_increment;

--  Change meta_coord.coord_system_id to unsigned
alter table meta_coord modify `coord_system_id` int(10) unsigned NOT NULL;


--  Change probe_design.coord_system_id to unsigned
alter table probe_design modify `coord_system_id` int(10) unsigned NOT NULL;




INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_b.sql|FK_consistency_tweaks');


