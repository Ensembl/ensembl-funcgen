-- Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
@header patch_91_92_j.sql - fastqc table
@desc   fastqc table
*/

drop table if exists phantom_peak;

create table phantom_peak (
  phantom_peak_id int(17) unsigned not null auto_increment,
  analysis_id  smallint(5) unsigned default null,
  alignment_id int(15) unsigned not null,
  num_reads    int(12) unsigned not null,
  est_frag_len        double default null,
  est_frag_len_2      double default null,
  est_frag_len_3      double default null,
  corr_est_frag_len   double default null,
  corr_est_frag_len_2 double default null,
  corr_est_frag_len_3 double default null,
  phantom_peak        int(17) unsigned not null,
  corr_phantom_peak   double  default null,
  argmin_corr      int(14) default null,
  min_corr         double  default null,
  nsc              double  default null,
  rsc              double  default null,
  quality_tag      int(14) default null,
  primary key (phantom_peak_id)
);

insert into phantom_peak (
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

ALTER TABLE phantom_peak ADD CONSTRAINT alignment_id_unique UNIQUE (alignment_id);

alter table phantom_peak add column run_failed boolean default false;
alter table phantom_peak add column error_message text;

drop table if exists alignment_qc_phantom_peak;
drop table if exists result_set_qc_phantom_peak;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_91_92_j.sql|phantom peak table');
