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
@header patch_90_91_zb.sql - 
@desc   
*/

rename table feature_set_qc_prop_reads_in_peaks to peak_calling_qc_prop_reads_in_peaks;

alter table peak_calling_qc_prop_reads_in_peaks change column feature_set_qc_prop_reads_in_peaks_id peak_calling_qc_prop_reads_in_peaks_id   int(10) unsigned NOT NULL;
alter table peak_calling_qc_prop_reads_in_peaks change column feature_set_id peak_calling_id int(10) unsigned NOT NULL;

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_90_91_zb.sql|');
