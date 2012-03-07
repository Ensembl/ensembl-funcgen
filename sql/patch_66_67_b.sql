/** 
@header patch_66_67_b.sql - regulatory_attribute.attribute_feature_idx
@desc   Add new index to facilitate RegulatoryFeatureAdaptor::fetch_all_by_attribute_feature
*/


-- High cardinality field first as we currently never want to query just on feature_table
ALTER table regulatory_attribute ADD KEY attribute_feature_idx (`attribute_feature_id`, `attribute_feature_table`);
analyze table regulatory_attribute;
optimize table regulatory_attribute;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_66_67_b.sql|regulatory_attribute.attribute_feature_idx');


