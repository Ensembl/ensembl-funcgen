-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

# patch_52_53_i.sql
#
# title: probe PRIMARY KEY
#
# description:
# Alter primary key to enable probes of same name but different array_chip i.e. Illumina

ALTER table probe add UNIQUE KEY `tmp_primary` (`probe_id`,`name`, `array_chip_id`);


ALTER table probe drop primary key;

ALTER table probe add PRIMARY KEY (`probe_id`,`name`, `array_chip_id`);
ALTER table probe drop key `tmp_primary`;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_i.sql|probe PRIMARY KEY');


