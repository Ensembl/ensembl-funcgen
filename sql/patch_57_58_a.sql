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

# patch_57_58_a.sql
#
# title: update schema version
#
# description:
# Update schema_version in meta table to 58

#We need to insert this for new species!
#Unique key will not capture duplicate NULL species_id patches
#So we need to handle this

#Is there a HC for multiple schema_version entries


#SET @schema_version_count=(SELECT COUNT(*) FROM meta WHERE meta_key='schema_version'); 
#This flow control only allows return values to be constants, not selects!
#Would need stored procedure/program/function with IF statement 
#IF (@schema_version_count > 1, select "Found more than 1 meta entry for schem_version. Please resolve", select NULL);
#IF (@schema_version_count = 1, select "updating", select NULL);	
#Or can we just delete and re-insert, living with the meta_id increment 

UPDATE meta SET meta_value='58' WHERE meta_key='schema_version';

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_57_58_a.sql|schema_version');


