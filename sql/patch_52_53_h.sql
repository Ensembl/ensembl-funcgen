# patch_52_53_c.sql
#
# title: feature_set.display_name
#
# description:
# Add feature_set.display_name

ALTER table feature_set add column `display_label` varchar(80) default NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_h.sql|feature_set.display_label');


