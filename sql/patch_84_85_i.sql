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
@header patch_84_85_i.sql - Normalise regulatory feature table: Create a non redundant version of the regulatory features
@desc   Also removes column "projected".
*/

CREATE TABLE regulatory_feature_nr (
  regulatory_feature_id int(10) unsigned NOT NULL auto_increment,
  feature_type_id int(10) unsigned default NULL,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_strand tinyint(1) NOT NULL,
  seq_region_start int(10) unsigned NOT NULL,
  seq_region_end int(10) unsigned NOT NULL,
  stable_id  varchar(18) DEFAULT NULL,
  bound_start_length mediumint(3) unsigned NOT NULL,
  bound_end_length mediumint(3) unsigned NOT NULL,
  epigenome_count smallint(6),
  PRIMARY KEY  (regulatory_feature_id),
  KEY feature_type_idx (feature_type_id),
  KEY stable_id_idx (stable_id)
) ENGINE=MyISAM;

insert into regulatory_feature_nr (
  seq_region_id,
  seq_region_start,
  seq_region_end,
  seq_region_strand,
  feature_type_id,
  stable_id,
  bound_start_length,
  bound_end_length,
  epigenome_count
) select 
  seq_region_id,
  seq_region_start,
  seq_region_end,
  seq_region_strand,
  feature_type_id,
  stable_id,
  bound_start_length,
  bound_end_length,
  epigenome_count
from regulatory_feature 
group by 
  seq_region_id,
  seq_region_start,
  seq_region_end,
  seq_region_strand,
  feature_type_id,
  stable_id,
  bound_start_length,
  bound_end_length,
  epigenome_count
;

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_i.sql|Normalise regulatory feature table: Create a non redundant version of the regulatory features.');


