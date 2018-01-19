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


/**
@header patch_75_76_g.sql - miRNA_target_feature
@desc Table definition  mirna_target_feature

*/

CREATE TABLE `mirna_target_feature` (
    `mirna_target_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
    `feature_set_id` int(10) unsigned NOT NULL,
    `feature_type_id` int(10) unsigned DEFAULT NULL,
    `accession` varchar(60) DEFAULT NULL,
    `display_label` varchar(60) DEFAULT NULL,
    `evidence` varchar(60) DEFAULT NULL,
    `interdb_stable_id` int(10) unsigned DEFAULT NULL,
    `method` varchar(60) DEFAULT NULL,
    `seq_region_id` int(10) unsigned NOT NULL,
    `seq_region_start` int(10) unsigned NOT NULL,
    `seq_region_end` int(10) unsigned NOT NULL,
    `seq_region_strand` tinyint(1) NOT NULL,
    `supporting_information` varchar(100) DEFAULT NULL,
    PRIMARY KEY (`mirna_target_feature_id`),
    UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
    KEY `feature_type_idx` (`feature_type_id`),
    KEY `feature_set_idx` (`feature_set_id`),
    KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_g.sql|mirna_target_feature');
