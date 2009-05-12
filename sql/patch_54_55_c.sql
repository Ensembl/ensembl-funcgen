# patch_54_55_c.sql
#
# title: coord_system.is_current_version_tidy
#
# description:
# Add an is_current field to the coord_system table to enable mart to filter the coord_systems
# Default is True as we assume newly added coord_system entries will always be current
# True clashes should be cleared up by update DB for release, or warned by the CoordSystem adaptor?
# Also clean up versions in non-assembled levels and redefine keys for multi species DBs



alter table coord_system add column is_current boolean default True;

alter table coord_system modify `version` varchar(255) DEFAULT NULL;

update coord_system set version=NULL where name not in ('scaffold', 'chromosome');


INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_c.sql|coord_system.is_current_version_null');

