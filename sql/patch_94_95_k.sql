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
@header patch_94_95_k.sql - peak calling statistic table
@desc   peak calling statistic table
*/

drop table if exists peak_calling_statistic_old;
create table peak_calling_statistic_old as select * from peak_calling_statistic;

drop table if exists peak_calling_statistic;
CREATE TABLE peak_calling_statistic (
  peak_calling_statistic_id int(28) unsigned NOT NULL AUTO_INCREMENT,
  peak_calling_id           int(18) unsigned DEFAULT NULL,
  epigenome_id              int(15) unsigned DEFAULT NULL,
  feature_type_id           int(18) unsigned DEFAULT NULL,
  statistic                 varchar(255) NOT NULL,
  value                     float(16) unsigned DEFAULT NULL,
  PRIMARY KEY (peak_calling_statistic_id)
);

insert into peak_calling_statistic(peak_calling_statistic_id, peak_calling_id, epigenome_id, feature_type_id, statistic, value)
select peak_calling_statistic_id, peak_calling_id, epigenome_id, feature_type_id, 'total_length', total_length from peak_calling_statistic_old join peak_calling using (peak_calling_id);

insert into peak_calling_statistic(peak_calling_id, epigenome_id, feature_type_id, statistic, value)
select  peak_calling_id, epigenome_id, feature_type_id, 'num_peaks', num_peaks from peak_calling_statistic_old join peak_calling using (peak_calling_id);

insert into peak_calling_statistic(peak_calling_id, epigenome_id, feature_type_id, statistic, value)
select peak_calling_id, epigenome_id, feature_type_id, 'average_length', average_length from peak_calling_statistic_old join peak_calling using (peak_calling_id);

drop table if exists peak_calling_statistic_old;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_94_95_k.sql|peak calling statistic table');
