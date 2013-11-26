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

# patch_52_53_f.sql
#
# title: array.class 
#
# description:
# Add class field to array table to define class of array e.g. AFFY_UTR, AFFY_ST etc (This sould be used to direct parsing). Format for expression arrays would be TRANSCRIPT.

ALTER table array add column `class` varchar(20) default NULL;
ALTER table array add UNIQUE KEY `class_name_idx` (`class`, `name`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_f.sql|array.class');


