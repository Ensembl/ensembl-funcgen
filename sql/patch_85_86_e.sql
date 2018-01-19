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
@header patch_85_86_e.sql - Add QC tables
@desc   Add QC tables
*/

CREATE TABLE `feature_set_qc_prop_reads_in_peaks` (
  `feature_set_qc_prop_reads_in_peaks_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id` int(10) unsigned DEFAULT NULL,
  `feature_set_id` int(10) unsigned NOT NULL,
  `prop_reads_in_peaks` double DEFAULT NULL,
  `total_reads` int(10) DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  `bam_file` varchar(512) NOT NULL,
  PRIMARY KEY (`feature_set_qc_prop_reads_in_peaks_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `result_set_qc_phantom_peak` (
  `result_set_qc_phantom_peak_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id` int(10) unsigned DEFAULT NULL,
  `result_set_id` int(10) unsigned NOT NULL,
  `filename` varchar(512) NOT NULL,
  `numReads` int(10) unsigned NOT NULL,
  `estFragLen` double DEFAULT NULL,
  `estFragLen2` double DEFAULT NULL,
  `estFragLen3` double DEFAULT NULL,
  `corr_estFragLen` double DEFAULT NULL,
  `corr_estFragLen2` double DEFAULT NULL,
  `corr_estFragLen3` double DEFAULT NULL,
  `phantomPeak` int(10) unsigned NOT NULL,
  `corr_phantomPeak` double DEFAULT NULL,
  `argmin_corr` int(10) DEFAULT NULL,
  `min_corr` double DEFAULT NULL,
  `NSC` double DEFAULT NULL,
  `RSC` double DEFAULT NULL,
  `QualityTag` int(10) DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  PRIMARY KEY (`result_set_qc_phantom_peak_id`),
  KEY `filename_idx` (`filename`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `result_set_qc_flagstats` (
  `result_set_qc_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `result_set_id` int(10) unsigned DEFAULT NULL,
  `analysis_id` int(10) unsigned DEFAULT NULL,
  `category` varchar(100) NOT NULL,
  `qc_passed_reads` int(10) unsigned DEFAULT NULL,
  `qc_failed_reads` int(10) unsigned DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  `bam_file` varchar(512) NOT NULL,
  PRIMARY KEY (`result_set_qc_id`),
  UNIQUE KEY `name_exp_idx` (`result_set_qc_id`,`category`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `result_set_qc_chance` (
  `result_set_qc_chance_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `signal_result_set_id` int(10) DEFAULT NULL,
  `analysis_id` int(10) unsigned DEFAULT NULL,
  `p` double DEFAULT NULL,
  `q` double DEFAULT NULL,
  `divergence` double DEFAULT NULL,
  `z_score` double DEFAULT NULL,
  `percent_genome_enriched` double DEFAULT NULL,
  `input_scaling_factor` double DEFAULT NULL,
  `differential_percentage_enrichment` double DEFAULT NULL,
  `control_enrichment_stronger_than_chip_at_bin` double DEFAULT NULL,
  `first_nonzero_bin_at` double DEFAULT NULL,
  `pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome` double DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  PRIMARY KEY (`result_set_qc_chance_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_85_86_e.sql|Add QC tables');
