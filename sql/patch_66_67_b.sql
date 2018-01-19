-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

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


