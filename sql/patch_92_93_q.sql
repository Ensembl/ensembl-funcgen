-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
@header patch_92_93_q.sql - Create segmentation_state_emission table
@desc   Create segmentation_state_emission table
*/

drop table if exists segmentation_state_emission;

CREATE TABLE segmentation_state_emission (
  segmentation_state_emission_id int(27) unsigned NOT NULL AUTO_INCREMENT,
  state    int(7) NOT NULL,
  CTCF     double NOT NULL,
  DNase1   double NOT NULL,
  H3K27ac  double NOT NULL,
  H3K27me3 double NOT NULL,
  H3K36me3 double NOT NULL,
  H3K4me1  double NOT NULL,
  H3K4me2  double NOT NULL,
  H3K4me3  double NOT NULL,
  H3K9ac   double NOT NULL,
  H3K9me3  double NOT NULL,
  PRIMARY KEY (segmentation_state_emission_id)
);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_q.sql|Create segmentation_state_emission table');
