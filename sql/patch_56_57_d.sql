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

# patch_56_57_d.sql
#
# title: CompareSchema tidy up
#
# description:
# Re-apply some patches which have been identified by the CompareSchema 
# as absent.

#Some array tables had wrong index name
ALTER IGNORE table array drop KEY `vendor`;
ALTER IGNORE table array drop KEY `vendor_name_idx`;
ALTER IGNORE table array add UNIQUE KEY `vendor_name_idx` (`vendor`, `name`);

#Lots of very minor column definition changes
ALTER table experimental_set modify `name` varchar(100) NOT NULL;
alter table result_set modify `analysis_id` smallint(5) unsigned NOT NULL;
ALTER table unmapped_object modify `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType', 'Probe', 'ProbeSet', 'ProbeFeature') NOT NULL;
ALTER table xref modify `description` VARCHAR(255);
alter table probe_design modify `analysis_id` smallint(5) unsigned NOT NULL;
alter table regulatory_feature modify `bound_seq_region_start` int(10) unsigned NOT NULL AFTER stable_id;
alter table regulatory_feature modify `bound_seq_region_end` int(10) unsigned NOT NULL AFTER bound_seq_region_start;
alter table coord_system modify `is_current` boolean default True AFTER species_id;
alter table probe_set modify `name` varchar(100) NOT NULL;
ALTER table cell_type change column description description varchar(80);
ALTER table external_db modify `type` ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL') default NULL;

#From patch_50_51_b.sql
ALTER IGNORE TABLE meta DROP INDEX key_value;
ALTER IGNORE TABLE meta DROP INDEX meta_key_index;
ALTER IGNORE TABLE meta DROP INDEX meta_value_index;
ALTER IGNORE TABLE meta DROP INDEX species_key_value_idx;
ALTER IGNORE TABLE meta DROP INDEX species_value_idx;
ALTER TABLE meta ADD UNIQUE INDEX species_key_value_idx (species_id, meta_key, meta_value);
ALTER TABLE meta ADD INDEX species_value_idx (species_id, meta_value);

#Dros only?
ALTER IGNORE TABLE meta_coord DROP INDEX cs_table_name_idx;
ALTER IGNORE TABLE meta_coord DROP INDEX table_name;
ALTER TABLE meta_coord ADD UNIQUE KEY `table_name` (`table_name`,`coord_system_id`);
ALTER IGNORE table probe_set ADD INDEX `name` (`name`);


#Correct result_feature index to remove PRIMARY KEY
#as it does not match the partition key
ALTER IGNORE table result_feature ADD KEY `result_feature_idx` (`result_feature_id`);
#remove PRIMARY KEY if present
ALTER IGNORE table result_feature DROP PRIMARY KEY; 
ALTER IGNORE table result_feature DROP INDEX result_feature_id; 
#Now need to drop other key and re instate as this will have dropped to the bottom?
#Is this really an issue, as it seems CompareSchema does not complain about
#the order of keys 
##ALTER table result_feature DROP INDEX set_window_seq_region_idx; 
##ALTER table result_feature ADD INDEX `set_window_seq_region_idx` (`result_set_id`, `window_size`,`seq_region_id`,`seq_region_start`);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_d.sql|CompareSchema_tidyup');


