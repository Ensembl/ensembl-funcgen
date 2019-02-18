-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
@header patch_96_97_d - Update mirna_target_feature
@desc   Removed link to FeatureSet, the only use was the link to Analysis,
        which has been added. Removed interdb_stable_id as never used
        Added Gene stable ID to simply link to gene
*/

ALTER TABLE mirna_target_feature DROP COLUMN interdb_stable_id;
ALTER TABLE mirna_target_feature ADD COLUMN analysis_id smallint(10) unsigned;
UPDATE mirna_target_feature SET analysis_id =  (SELECT analysis_id FROM feature_set WHERE name = 'TarBase miRNA');
ALTER TABLE mirna_target_feature DROP COLUMN feature_set_id;
ALTER TABLE mirna_target_feature ADD COLUMN gene_stable_id VARCHAR(128);
UPDATE mirna_target_feature mtf
  JOIN object_xref ox ON ox.ensembl_id = mtf.mirna_target_feature_id AND
       ensembl_object_type = 'MirnaTargetFeature' JOIN xref USING (xref_id)
  SET gene_stable_id = dbprimary_acc;
DELETE object_xref, xref FROM object_xref JOIN xref USING (xref_id) WHERE ensembl_object_type='MirnaTargetFeature';
ALTER TABLE `mirna_target_feature` ADD UNIQUE KEY `unique_idx` (`accession`,`gene_stable_id`,`seq_region_start`,`seq_region_end`);

-- patch identifier

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_96_97_d.sql|Update mirna_target_feature');
