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
@header patch_69_70_b.sql - regulatory_feature.bound_start/end_length
@desc   Change the bound_seq_regiobn_start/end fields to lengths
*/

ALTER TABLE regulatory_feature ADD `bound_start_length` mediumint(3) unsigned NOT NULL;
ALTER TABLE regulatory_feature ADD `bound_end_length` mediumint(3) unsigned NOT NULL;

update regulatory_feature set bound_start_length=(seq_region_start - bound_seq_region_start);
update regulatory_feature set bound_end_length=(bound_seq_region_end - seq_region_end);

ALTER TABLE regulatory_feature DROP bound_seq_region_start;
ALTER TABLE regulatory_feature DROP bound_seq_region_end;

analyze table regulatory_feature;
optimize table regulatory_feature;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_69_70_b.sql|regulatory_feature.bound_start/end_length');

