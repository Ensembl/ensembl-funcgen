-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
@header patch_92_93_v.sql - Create probemapping statistic table
@desc   Create probemapping statistic table
*/

drop table if exists probe_mapping_statistic;

CREATE TABLE probe_mapping_statistic (
  probe_mapping_statistic_id int(29) unsigned NOT NULL AUTO_INCREMENT,
  array_id  int(11) unsigned DEFAULT NULL,
  statistic varchar(255) NOT NULL,
  value     double unsigned DEFAULT NULL,
  PRIMARY KEY (probe_mapping_statistic_id)
);

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'number_of_probes', count(distinct probe_id) from array join array_chip using (array_id) join probe using (array_chip_id) group by array.name order by array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'number_of_probe_sets', count(distinct probe_set_id) from array join array_chip using (array_id) join probe_set using (array_chip_id) group by array.name order by array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'number_of_probe_features', count(distinct probe_feature.probe_feature_id) from array join array_chip using (array_id) join probe using (array_chip_id) join probe_feature using (probe_id) group by array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'number_of_probe_features_from_transcripts', count(distinct probe_feature.probe_feature_id) from array join array_chip using (array_id) join probe using (array_chip_id) join probe_feature using (probe_id) join analysis using (analysis_id) where analysis.logic_name like "%transcript%" group by analysis.logic_name, array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array_chip.array_id, 'number_of_unmapped_probes', count(probe.probe_id) from array_chip join probe using (array_chip_id) left join probe_feature using (probe_id) where probe_feature.probe_id is null group by array_chip.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'average_probe_length', avg(length(sequence)) from array join array_chip using (array_id) join probe using (array_chip_id) join probe_seq using (probe_seq_id) group by array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'max_probe_length', max(length(sequence)) from array join array_chip using (array_id) join probe using (array_chip_id) join probe_seq using (probe_seq_id) group by array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'min_probe_length', min(length(sequence)) from array join array_chip using (array_id) join probe using (array_chip_id) join probe_seq using (probe_seq_id) group by array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'number_of_probes_mapped_to_transcripts', count(distinct probe_transcript.probe_id) from array join array_chip using (array_id) join probe using (array_chip_id) left join probe_transcript using (probe_id) group by array.array_id;

-- Too slow!
-- insert into probe_mapping_statistic (array_id, statistic, value) 
-- select array.array_id, 'number_of_probe_sets_mapped_to_transcripts', count(distinct probe_set_transcript.probe_set_id) from array join array_chip using (array_id) join probe_set using (array_chip_id) left join probe_set_transcript using (probe_set_id) group by array.array_id;

insert into probe_mapping_statistic (array_id, statistic, value) 
select array.array_id, 'number_of_probe_sets', count(distinct probe_set.probe_set_id) from array join array_chip using (array_id) left join probe_set using (array_chip_id) group by array.name order by array.array_id;


-- patch identifier
INSERT INTO `meta` (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_v.sql|Create probemapping statistic table');
