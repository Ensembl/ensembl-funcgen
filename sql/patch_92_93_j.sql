-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
@header patch_91_92_j.sql - Add table to store phantom peak outcomes
@desc   Add table to store phantom peak outcomes
*/

drop table if exists phantom_peak;

CREATE TABLE phantom_peak (
  phantom_peak_id int(17) unsigned NOT NULL AUTO_INCREMENT,
  analysis_id smallint(5) unsigned DEFAULT NULL,
  alignment_id int(15) unsigned NOT NULL,
  num_reads int(12) unsigned NOT NULL,
  est_frag_len double DEFAULT NULL,
  est_frag_len_2 double DEFAULT NULL,
  est_frag_len_3 double DEFAULT NULL,
  corr_est_frag_len double DEFAULT NULL,
  corr_est_frag_len_2 double DEFAULT NULL,
  corr_est_frag_len_3 double DEFAULT NULL,
  phantom_peak int(17) unsigned NOT NULL,
  corr_phantom_peak double DEFAULT NULL,
  argmin_corr int(14) DEFAULT NULL,
  min_corr double DEFAULT NULL,
  nsc double DEFAULT NULL,
  rsc double DEFAULT NULL,
  quality_tag int(14) DEFAULT NULL,
  run_failed tinyint(1) DEFAULT '0',
  error_message text,
  PRIMARY KEY (phantom_peak_id),
  UNIQUE KEY alignment_id_unique (alignment_id)
);

insert ignore into phantom_peak (
  analysis_id,
  alignment_id,
  num_reads,
  est_frag_len,
  est_frag_len_2,
  est_frag_len_3,
  corr_est_frag_len,
  corr_est_frag_len_2,
  corr_est_frag_len_3,
  phantom_peak,
  corr_phantom_peak,
  argmin_corr,
  min_corr,
  nsc,
  rsc,
  quality_tag
)
  select
    analysis_id,
    alignment_id,
    numReads,
    estFragLen,
    estFragLen2,
    estFragLen3,
    corr_estFragLen,
    corr_estFragLen2,
    corr_estFragLen3,
    phantomPeak,
    corr_phantomPeak,
    argmin_corr,
    min_corr,
    NSC,
    RSC,
    QualityTag
  from 
    alignment_qc_phantom_peak
;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_j.sql|phantom peak table');
