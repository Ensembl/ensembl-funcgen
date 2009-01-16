# patch_52_53_c.sql
#
# title: regulatory_feature.bound_seq_region_start/end
#
# description:
# Add bound_seq_region_start and end to regualtory feature table to 
# prevent having to calculate on the fly

ALTER table regulatory_feature add column `bound_seq_region_start` int(10) unsigned NOT NULL;
ALTER table regulatory_feature add column `bound_seq_region_end` int(10) unsigned NOT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_g.sql|regulatory_feature.bound_seq_region_start/end');


