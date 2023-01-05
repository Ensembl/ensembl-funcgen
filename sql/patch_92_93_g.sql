-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
@header patch_91_92_g.sql - New table idr
@desc   New table idr
*/

DROP TABLE IF EXISTS idr;
CREATE TABLE idr (
  idr_id int(9) unsigned NOT NULL AUTO_INCREMENT,
  experiment_id int(15) unsigned NOT NULL,
  max_peaks int(11) unsigned NOT NULL,
  type enum('on biological replicates','on technical replicates') NOT NULL,
  PRIMARY KEY (idr_id)
);

-- alter table idr add column failed_idr_pairs varchar(255) default null;
alter table idr add column failed_idr_pairs text default null;
alter table idr change column type type enum('on biological replicates','on technical replicates', 'no_idr') NOT NULL;
alter table idr change column max_peaks max_peaks int(11) unsigned default NULL;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_g.sql|New table idr');
