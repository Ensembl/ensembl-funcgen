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
@header patch_91_92_d.sql - Create table for chance quality check
@desc   Create table for chance quality check
*/

drop table if exists alignment_qc_chance;

CREATE TABLE if not exists `chance` (
  `chance_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `signal_alignment_id` int(10),
  `control_alignment_id` int(10),
  `analysis_id`        int(10) unsigned,
-- Not really that important
-- See slide 38 on 
-- http://www.ebi.ac.uk/seqdb/confluence/download/attachments/18483313/UCL_ChIPseq_Wilder.pptx?version=2&modificationDate=1442910347000&api=v2
-- dashed green line
--
  `p` double default NULL,
-- Not really that important
  `q` double default NULL,
--
-- This is the main statistic.
-- It is a scaled version of differential_percentage_enrichment. The reason 
-- is that the exact location is important and that is not reflected in 
-- differential_percentage_enrichment.
-- 
--
  `divergence` double default NULL,
--
-- Distance from the mean, if the distribution was standardised to a normal distribution
--
  `z_score` double default NULL,
--
-- Distance between dashed green line and 1
--
  `percent_genome_enriched` double default NULL,
--
-- A suggestion on how to scale the control to equal the background noise 
-- in the signal
--
  `input_scaling_factor` double default NULL,
--
-- It is the greates distance between the cumulative coverage lines of the 
-- control and the signal when plotted into a graph.
--
  `differential_percentage_enrichment` double default NULL,
--
-- Usually the two curves would meet at one. If there is an enrichment in 
-- the control, then this reports the bin number where this happens.
-- A bit visible in diagram d, slide 40 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053734/figure/F2/
--
  `control_enrichment_stronger_than_chip_at_bin`double default NULL,
-- 
-- After sorting the bins from the signal, this is the rank of the first non zero bin.
--
  `first_nonzero_bin_at`double default NULL,
--
-- Proportion of control reads in the highest 1 percent of the bins. The expected value would be 0.01, but only
-- greater deviations from that are reported.
--
  `pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome`double default NULL,
  `path` varchar(512) NOT NULL,
  PRIMARY KEY (`chance_id`)
);

ALTER TABLE chance ADD CONSTRAINT signal_control_alignment_unique UNIQUE (signal_alignment_id, control_alignment_id);

alter table chance add column run_failed boolean default false;
alter table chance add column error_message text;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_d.sql|Create table for chance quality check');
