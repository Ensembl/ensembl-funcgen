-- Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
@header patch_91_92_h.sql - execution_plan table
@desc   execution_plan table
*/

DROP TABLE IF EXISTS execution_plan;

CREATE TABLE execution_plan (
  execution_plan_id int(18) unsigned NOT NULL AUTO_INCREMENT,
  time BIGINT default null,
  experiment_id int(16) unsigned NOT NULL,
  execution_plan text,
  PRIMARY KEY (execution_plan_id)
);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_h.sql|execution_plan table');
