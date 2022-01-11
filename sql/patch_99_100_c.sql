-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
@header patch_99_100_c.sql - remove default values from text data type
@desc   The TEXT data type can not be assigned a default value
*/

ALTER TABLE regulatory_build MODIFY COLUMN description TEXT;
ALTER TABLE epigenome MODIFY COLUMN description MEDIUMTEXT;
ALTER TABLE epigenome MODIFY COLUMN search_terms MEDIUMTEXT;
ALTER TABLE epigenome MODIFY COLUMN full_name MEDIUMTEXT;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_99_100_c.sql|remove default values from text data type');
