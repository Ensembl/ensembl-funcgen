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


ALTER TABLE coord_system ADD `core_coord_system_id` int(10) NOT NULL;
drop index rank on coord_system;
drop index name on coord_system;
ALTER TABLE coord_system  ADD KEY `name_version_idx` (`name`, `version`);
ALTER TABLE coord_system change schema_build schema_build varchar(8) default NULL;

select "You now need to manually update the core_coord_system_ids and resolve any name version conflicts before continuing";
exit;


ALTER TABLE coord_system DROP PRIMARY KEY, ADD PRIMARY KEY(`coord_system_id`, `core_coord_system_id`, `schema_build`);

select "All name version pairs should correspond to the same(nr) coord_system_id, now update your feature and meta tables if required";
exit;

-- array.type and new key
alter table array add UNIQUE KEY   (`vendor`, `name`)
alter table array add `type` varchar(20) default NULL;
select "You need to manually update your array.type as OLIGO or PCR";



-- probe array_chip_idx
alter table probe add KEY `array_chip_idx` (`array_chip_id`);


-- probe_feature.cigar_line
-- alter table probe_feature add `cigar_line` text;
