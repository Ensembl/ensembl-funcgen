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

# patch_53_54_c.sql
#
# title: Add multispecies support
#
# description:
# Modifies coord_system and seq_region to ensure that the primary keys and structure allows for multispecies 



# Auto-inc fields are indexed before primary key removal to avoid an ERROR 1075 (42000).

ALTER TABLE coord_system ADD column `species_id` int(10) DEFAULT 1;
ALTER TABLE coord_system ADD index `coord_species_idx` (`species_id`);
ALTER TABLE coord_system ADD index `coord_system_id_idx` (`coord_system_id`);
ALTER TABLE coord_system DROP primary key;
ALTER TABLE coord_system ADD primary key (`name`,`version`,`schema_build`,`species_id`);

ALTER TABLE seq_region ADD index `seq_region_id_idx` (`seq_region_id`);
ALTER TABLE seq_region DROP PRIMARY KEY;
ALTER TABLE seq_region ADD PRIMARY KEY  (`name`, `schema_build`, `coord_system_id`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_53_54_c.sql|multispecies');
