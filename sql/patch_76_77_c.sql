-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
@header patch_76_77_c.sql - Correct FeatureType mirna so_name and accession
@desc   so_name and accessions have been swapped in some cases, this patch fixes these
*/

UPDATE feature_type SET so_accession = 'SO:0000934'        WHERE so_accession = 'miRNA_target_site';
UPDATE feature_type SET so_name      = 'miRNA_target_site' WHERE so_name      = 'SO:0000934';

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_76_77_c.sql|Correct mirna so_name and accession in feature_type');



