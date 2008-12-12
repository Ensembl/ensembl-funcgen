# patch_52_53_d.sql
#
# title: unmapped_object ensembl_object_type CoordSystem
#
# description:
# Modify unmapped_object_type.ensembl_object_type to include CoordSystem

ALTER table unmapped_object modify `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType', 'Probe', 'CoordSystem') NOT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_d.sql|ensembl_object_type CoordSystem');


