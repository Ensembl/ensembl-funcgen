-- Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


