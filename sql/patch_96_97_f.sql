-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
@header patch_96_97_f.sql - Add search_terms and full_name columns to epigenome table and rename display_label column to short_name
<<<<<<< HEAD
@desc 	Add search_terms and full_name columns to epigenome table, rename display_label column to short_name and change description to TEXT
=======
@desc 	Add search_terms and full_name columns to epigenome table and rename display_label column to short_name
>>>>>>> 122d95d901ba665f2dc9e9d6ae180b54567a7ba2
*/

ALTER TABLE epigenome ADD COLUMN search_terms MEDIUMTEXT DEFAULT NULL;
ALTER TABLE epigenome ADD COLUMN full_name MEDIUMTEXT DEFAULT NULL;
ALTER TABLE epigenome DROP INDEX display_label_idx;
ALTER TABLE epigenome CHANGE `display_label` `short_name` varchar(120) NOT NULL;
ALTER TABLE epigenome ADD CONSTRAINT short_name_idx UNIQUE (short_name);
ALTER TABLE epigenome MODIFY description MEDIUMTEXT DEFAULT NULL;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_96_97_f.sql|Add search_terms and full_name columns to epigenome table, rename display_label column to short_name and change description to TEXT');

