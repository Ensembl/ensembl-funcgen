# patch_52_53_c.sql
#
# title: ensembl_object_type Probe 
#
# description:
# Modify object_xref.ensembl_object_type to include Probe

ALTER table object_xref modify `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType', 'ProbeSet', 'Probe') NOT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_c.sql|ensembl_object_type ProbeSet Probe');


