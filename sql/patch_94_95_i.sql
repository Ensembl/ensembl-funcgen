-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
@header patch_94_95_i.sql - segmentation_statistic table
@desc   segmentation_statistic table
*/

drop table if exists segmentation_statistic;
CREATE TABLE segmentation_statistic (
  segmentation_statistic_id int(30) unsigned NOT NULL AUTO_INCREMENT,
  segmentation_id           int(18) unsigned default null,
  state                     int(8)  unsigned default null,
  epigenome_id              int(22) unsigned default NULL,
  label                     varchar(255) default NULL,
  statistic                 varchar(255) NOT NULL,
  value                     float(11) unsigned DEFAULT null,
  PRIMARY KEY (segmentation_statistic_id),
  UNIQUE KEY stats_uniq (statistic,segmentation_id,epigenome_id,label)
);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_94_95_i.sql|segmentation_statistic table');
