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

# patch_56_57_j.sql
#
# title: probe_feature.cigar_line_format
#
# description:
# Update cigar_line format to extended format agreed by SAM-Tools group


#Seq & Alignment match
update probe_feature set cigar_line=replace(cigar_line, 'M', '=');
update probe_feature set cigar_line=replace(cigar_line, 'V', '=');

#Alignment match & Seq mismatch
update probe_feature set cigar_line=replace(cigar_line, 'm', 'X');

#Soft clipping - Used for unknown sequence overhanging the end of a transcript
update probe_feature set cigar_line=replace(cigar_line, 'U', 'S');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_j.sql|probe_feature.cigar_line_format');


 
