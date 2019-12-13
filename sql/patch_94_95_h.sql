-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
@header patch_94_95_h.sql - segmentation_cell_tables
@desc   segmentation_cell_tables
*/

drop table if exists segmentation_cell_tables;
CREATE TABLE if not exists segmentation_cell_tables (
  superclass           varchar(255) NOT NULL,
  class                varchar(255) NOT NULL,
  segmentation_id      int(18) unsigned NOT NULL,
  epigenome_id         int(16) unsigned NOT NULL,
  feature_type_id      int(18) unsigned NOT NULL,
  signal_alignment_id  int(23) unsigned NOT NULL,
  control_alignment_id int(23) unsigned NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_94_95_h.sql|segmentation_cell_tables');
