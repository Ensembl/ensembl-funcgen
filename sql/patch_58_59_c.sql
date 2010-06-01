# patch_58_59_c.sql
#
# title: add annotated_feature.summit
#
# description:
# Add summit field to define peak summit of an annotated_feature


ALTER table annotated_feature ADD `summit` int(10) unsigned default NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_58_59_c.sql|annotated_feature.summit');


