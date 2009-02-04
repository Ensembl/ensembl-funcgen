# patch_52_53_j.sql
#
# title: add ProbeFeature to xref schema ensembl_object_type enums
#
# description:
# Modify identity_xref changing 

ALTER table object_xref modify `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType','ProbeSet','Probe', 'ProbeFeature') NOT NULL;
ALTER table unmapped_object modify `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType','ProbeSet','Probe', 'ProbeFeature') NOT NULL; 



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_j.sql|ensembl_object_type ProbeFeature');


