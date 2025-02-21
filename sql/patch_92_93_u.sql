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
@header patch_92_93_u.sql - Create probemapping meta table
@desc   Create probemapping meta table
*/

drop table if exists probe_mapping;

CREATE TABLE probe_mapping (
  probe_mapping_id    int(22) unsigned NOT NULL AUTO_INCREMENT,
  assembly            varchar(255) DEFAULT NULL,
  gene_build_version  varchar(255) DEFAULT NULL,
  five_prime_utr      int(22) unsigned DEFAULT NULL,
  three_prime_utr     int(22) unsigned DEFAULT NULL,
  sample_probe_id     int(22) unsigned DEFAULT NULL,
  sample_probe_set_id int(22) unsigned DEFAULT NULL,
  release_version     varchar(255) DEFAULT NULL,
  release_date        varchar(255) DEFAULT NULL,
  PRIMARY KEY (probe_mapping_id)
);

-- patch identifier
INSERT INTO `meta` (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_u.sql|Create probemapping meta table');
