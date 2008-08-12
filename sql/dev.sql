-- List features above 2kb.



--SELECT (t.seq_region_end - t.seq_region_start + 1 ) as max, t.annotated_feature_id, t.feature_set_id FROM annotated_feature t where max >2000 order by max desc limit 100;

-- Can't use alias in where/order

-- Ignore eQTLs, these should really be in a different table if they are likely to be much larger than other features
-- This will speed up the range queries

 SELECT (t.seq_region_end - t.seq_region_start + 1 ) as max, t.annotated_feature_id, t.feature_set_id FROM annotated_feature t where (t.seq_region_end - t.seq_region_start + 1 ) >2000 and t.feature_set_id!=70 order by (t.seq_region_end - t.seq_region_start + 1 ) desc;

-- get feature_sets effected
-- SELECT distinct(fs.name), t.annotated_feature_id, t.feature_set_id FROM annotated_feature t where (t.seq_region_end - t.seq_region_start + 1 ) >2000 order by (t.seq_region_end - t.seq_region_start + 1 ) desc;



-- Cannot do this for probe feature as we have the sanger PCR array, need to do a not in for the chip_cannel_ids.


-----------------------------------

--- associated_feature_types
--- should really change associated_feature_type to use a number for feature_table?
--- This is only true if we are actually matching the field on a query?
--- This will be thes case when we are retrieving using a Feature
--- But not when we are retrieving using a FeatureType
--- Main usage id by Feature, so number would be better
--- Not the case for RegFeat as we use the reg_feat_id as the key, not the table_id/name
--- We need an extra index on feature_type_id.
--- Don't want to obfuscate table_names in DB, maybe we could have them in meta?
--- And get the SetFeature adaptor to load them when instantiated
--- We would lose some speed on converting numbers to text when retrieving reg_attrs
--- Or pulling features back by their group feature_types

--- Is this always going to represent groups, in the sense of a motif group?
--- Or can this be used generically as simply a group of feature_types?
--- This is simply associating multiple feature types with the same feature
--- Maybe feature_type_set, no confusion with other sets.
--- feature_attribute? Does not have feature_type focus and maybe conflict with the idea of a core attribute.
--- Which 

--- groups should not really be added as FeatureTypes as they are a little more abstract.
--- Less of a region on the genome and more of a classification of a group of feature_types.


--- regulatory_attribute should be renamed regulatory_supporting_feature?
--- Then we could have annotated_supporting_feature..if required, or do we just consider all combined
--- features a regulatory_feature, and make the distinction with the feature set?

--- in core terms attrs can be a count of somthing e.g has a value, can be 1 or another metric.

--- So we're really talking about classes here, which may overlap, kinda like a DAG.
--- Then we have ambiguity with feature_type.class

--- Just implement this for now and think about a more generic solution?
--- methods would be get_FeatureTypes? get_attribute_FeatureTypes?
--- get_group_FeatureTypes get_associated_FeatureTypes?
