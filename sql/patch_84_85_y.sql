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
@header patch_84_85_y.sql - Table for storing epigenomes used in the regulatory build
@desc Table for storing epigenomes used in the regulatory build
*/

drop table if exists regulatory_build_epigenome;

create table regulatory_build_epigenome (
  regulatory_build_epigenome_id int(10) unsigned NOT NULL auto_increment,
  regulatory_build_id int(10) unsigned NOT NULL,
  epigenome_id int(10) unsigned NOT NULL,
  PRIMARY KEY  (regulatory_build_epigenome_id)
);

insert into regulatory_build_epigenome(
regulatory_build_id,
epigenome_id
)
select 
 distinct regulatory_build.regulatory_build_id, epigenome.epigenome_id
from 
 epigenome 
 join regulatory_activity using (epigenome_id) 
 join regulatory_feature using (regulatory_feature_id) 
 join regulatory_build using (regulatory_build_id);

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_y.sql|Table for storing epigenomes used in the regulatory build');
