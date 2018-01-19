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
@header patch_84_85_v.sql - Move meta entries regarding regulatory build to the regulatory_build table
@desc Move meta entries regarding regulatory build to the regulatory_build table
*/

update regulatory_build, meta m1, meta m2, meta m3 
set 
  initial_release_date=m1.meta_value,
  version=m2.meta_value,
  last_annotation_update=m3.meta_value
where 
  m1.meta_key="regbuild.initial_release_date" 
  and m2.meta_key="regbuild.version" 
  and m3.meta_key="regbuild.last_annotation_update";

delete from meta where meta_key in ("regbuild.initial_release_date", "regbuild.version", "regbuild.last_annotation_update");

insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_v.sql|Move meta entries regarding regulatory build to the regulatory_build table');
