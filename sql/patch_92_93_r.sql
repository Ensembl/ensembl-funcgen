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
@header patch_92_93_r.sql - Create regulatory_build_statistics table
@desc   Create regulatory_build_statistics table
*/

drop table if exists regulatory_build_statistic;

CREATE TABLE regulatory_build_statistic (
  regulatory_build_statistics_id int(30) unsigned NOT NULL AUTO_INCREMENT,
  regulatory_build_id int(22) unsigned DEFAULT NULL,
  statistic           varchar(255) DEFAULT NULL,
  value               double unsigned DEFAULT NULL,
  PRIMARY KEY (regulatory_build_statistics_id),
  UNIQUE KEY stats_uniq (statistic,regulatory_build_id)
);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_r.sql|Create regulatory_build_statistics table');
