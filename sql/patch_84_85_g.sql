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
@header patch_84_85_g.sql - update dbentry related tables.
@desc   Updates to catch up with developments from the core schema and allow xrefs to be stored for epigenomes.
*/

alter table external_db modify external_db_id INTEGER UNSIGNED NOT NULL AUTO_INCREMENT;
alter table external_db modify type ENUM(
  'ARRAY',
  'ALT_TRANS',
  'ALT_GENE',
  'MISC',
  'LIT',
  'PRIMARY_DB_SYNONYM', 
  'ENSEMBL'
);

alter table unmapped_object modify external_db_id INTEGER UNSIGNED;
alter table object_xref     modify ensembl_object_type ENUM (
  'Epigenome', 
  'Experiment', 
  'RegulatoryFeature', 
  'ExternalFeature', 
  'AnnotatedFeature', 
  'FeatureType', 
  'MirnaTargetFeature',
  'ProbeSet',
  'Probe',
  'ProbeFeature'
) not NULL;
alter table xref modify external_db_id INTEGER UNSIGNED;
alter table xref modify dbprimary_acc VARCHAR(512) NOT NULL;
alter table xref modify display_label VARCHAR(512) NOT NULL;
alter table xref modify version VARCHAR(10) DEFAULT NULL;
alter table xref modify description TEXT;
alter table xref modify info_type ENUM (
  'NONE', 
  'PROJECTION', 
  'MISC', 
  'DEPENDENT',
  'DIRECT', 
  'SEQUENCE_MATCH',
  'INFERRED_PAIR', 
  'PROBE',
  'UNMAPPED', 
  'COORDINATE_OVERLAP', 
  'CHECKSUM' 
) DEFAULT 'NONE' NOT NULL;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_84_85_g.sql|update dbentry related tables');
