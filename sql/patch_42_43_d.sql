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

-- update coord_system_ids and meta_coord and meta table
-- This will disappear when we implement version to assembly mapping 
-- Tidy up analysis table and ids;

-- update meta table
delete from meta where meta_key="schema.version";
insert into meta values('', 'schema_version', '43');

-- update coord_system_table dependent on species
select "Need to manually update the schema_build in coord_system to match the species schema_build";
--update coord_system set schema_build='43_36e' where schema_build='42_36d';
--update coord_system set schema_build='43_36d' where schema_build='42_36c';



-- Tidy up analysis tables/ids
update predicted_feature pf, analysis a set pf.analysis_id=a.analysis_id where a.logic_name='Nessie';
delete a, ad from analysis a, analysis_description ad where a.logic_name='TilingHMM' and a.analysis_id=ad.analysis_id ;
