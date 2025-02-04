-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
@header patch_104_105_b.sql - New indices
@desc   These should improve track loading on the genome browser
*/

CREATE INDEX peak_track ON peak (peak_calling_id, seq_region_id, seq_region_start, seq_region_end, score, seq_region_strand, summit, peak_id);
CREATE INDEX speedup ON motif_feature_regulatory_feature (regulatory_feature_id ASC, has_matching_Peak, motif_feature_id);

-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_104_105_b.sql|New indices');
