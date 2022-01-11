-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
@header patch_93_94_d.sql - Adds table segmentation_cell_table_without_ctcf
@desc   Adds table segmentation_cell_table_without_ctcf
*/

DROP TABLE IF EXISTS segmentation_cell_table_without_ctcf;

CREATE TABLE segmentation_cell_table_without_ctcf (
  epigenome        varchar(120) DEFAULT NULL,
  feature_type     varchar(40)  NOT NULL,
  signal_bam_path  varchar(255) NOT NULL,
  control_bam_path varchar(255) NOT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_93_94_d.sql|Adds table segmentation_cell_table_without_ctcf');
