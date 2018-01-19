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

# patch_59_60_f.sql
#
# title: bm.frequencies_pf_index_mods
#
# description:
# Make some minor modifications to the binding_matric frequences column
# and probe_feature indexes to increase performance

alter table binding_matrix modify `frequencies` varchar(1000) NOT NULL;

alter table probe_feature add index `seq_region_probe_probe_feature_idx` (`seq_region_id`,`seq_region_start`, `seq_region_end`, `probe_id`, `probe_feature_id`); 
alter table probe_feature drop index seq_region_idx; 


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_59_60_f.sql|bm.frequencies_pf.index_mods');


