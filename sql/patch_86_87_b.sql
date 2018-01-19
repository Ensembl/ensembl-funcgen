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
@header patch_86_87_b.sql - Change data type of certain columns to facilitate foreing key constraints
@desc   Change data type of certain columns to facilitate foreing key constraints
*/

ALTER TABLE regulatory_build MODIFY analysis_id smallint(5) unsigned NOT NULL;
ALTER TABLE feature_set_qc_prop_reads_in_peaks MODIFY analysis_id smallint(5) unsigned DEFAULT NULL;
ALTER TABLE result_set_qc_flagstats MODIFY analysis_id smallint(5) unsigned DEFAULT NULL;
ALTER TABLE result_set_qc_chance MODIFY analysis_id smallint(5) unsigned DEFAULT NULL;
ALTER TABLE result_set_qc_chance MODIFY signal_result_set_id int(10) unsigned DEFAULT NULL;
ALTER TABLE result_set_qc_phantom_peak MODIFY analysis_id smallint(5) unsigned DEFAULT NULL;
ALTER TABLE segmentation_file MODIFY regulatory_build_id int(4) unsigned DEFAULT NULL;


-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_86_87_b.sql|Change data type of certain columns to facilitate foreing key constraints');