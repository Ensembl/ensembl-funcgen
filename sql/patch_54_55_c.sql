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

# patch_54_55_c.sql
#
# title: coord_system.is_current_version_tidy
#
# description:
# Add an is_current field to the coord_system table to enable mart to filter the coord_systems
# Default is True as we assume newly added coord_system entries will always be current
# True clashes should be cleared up by update DB for release, or warned by the CoordSystem adaptor?


alter table coord_system add column is_current boolean default True;

alter table coord_system modify `version` varchar(255) DEFAULT NULL;


INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_c.sql|coord_system.is_current_version_null');

