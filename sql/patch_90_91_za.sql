-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
@header patch_90_91_za.sql - Move peak_callings from feature_set to peak_calling
@desc   Move peak_callings from feature_set to peak_calling
*/

DROP TABLE IF EXISTS `peak_calling`;

CREATE TABLE `peak_calling` (
  `peak_calling_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(300) NOT NULL,
  `display_label` varchar(300) NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  `analysis_id`     smallint(5) unsigned NOT NULL,
  `alignment_id`    int(10) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `experiment_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`peak_calling_id`),
  UNIQUE KEY `peak_calling_id_idx` (`peak_calling_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into peak_calling (
  peak_calling_id,
  name,
  display_label,
  feature_type_id,
  analysis_id,
  alignment_id,
  epigenome_id,
  experiment_id
) select 
  feature_set_id, 
  feature_set.name,
  feature_set.display_label,
  feature_type_id, 
  analysis_id, 
  supporting_set_id,
  epigenome_id,
  experiment_id
from 
  feature_set 
  join data_set using (feature_set_id) 
  join supporting_set using (data_set_id) 
where 
  feature_set.type="annotated";

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_za.sql|Move peak_callings from feature_set to peak_calling');
