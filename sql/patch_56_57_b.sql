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

# patch_56_57_b.sql
#
# title: unmapped_object object_type key
#
# description:
# Add a ensembl_object_type/id key to the unmapped_object table

#ALTER TABLE unmapped_object DROP KEY `object_type_idx`;

ALTER TABLE unmapped_object ADD KEY `object_type_idx` (`ensembl_id`, `ensembl_object_type`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_b.sql|uo.object_type_id_key');


