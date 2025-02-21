-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
@header patch_92_93_o.sql - Create peak_calling_statistic table
@desc   Create peak_calling_statistic table
*/

drop table if exists peak_calling_statistic;

CREATE TABLE peak_calling_statistic (
  peak_calling_statistic_id int(28) unsigned  NOT NULL AUTO_INCREMENT,
  peak_calling_id int(18)    NOT NULL DEFAULT '0',
  total_length    int(15)    DEFAULT NULL,
  num_peaks       bigint(14) DEFAULT NULL,
  average_length  int(17)    DEFAULT NULL,
  PRIMARY KEY (`peak_calling_statistic_id`)
);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_o.sql|Create peak_calling_statistic table');
