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

# patch_54_55_a.sql
#
# title: FK consistency tweaks
#
# description:
# Some general schema tweaks to bring var types in line with core and internally consistent
# To allow valid foreign keys when using innodb

-- Change analysis_id to smallint(5)
alter table analysis modify `analysis_id` smallint(5) unsigned NOT NULL auto_increment;
alter table probe_feature modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table feature_set modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table result_set modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table probe_design modify `analysis_id` smallint(5) unsigned NOT NULL;

-- Change s.status_name_id to unsigned
alter table status modify `status_name_id` int(10) unsigned NOT NULL;

--  Change c.coord_system_id to unsigned
alter table coord_system modify `coord_system_id` int(10) unsigned NOT NULL auto_increment;

--  Change meta_coord.coord_system_id to unsigned
alter table meta_coord modify `coord_system_id` int(10) unsigned NOT NULL;


--  Change probe_design.coord_system_id to unsigned
alter table probe_design modify `coord_system_id` int(10) unsigned NOT NULL;




INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_b.sql|FK_consistency_tweaks');


