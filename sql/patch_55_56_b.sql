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

# patch_55_56_b.sql
#
# title: probe_feature.cigar_line format
#
# description:
# Change cigar_line symbols to reflect those agreed in the extended SAM cigar line format

#Dont need COLLATE latin1_general_cs as replace performs case sensitive match

update probe_feature set cigar_line=replace(cigar_line, 'U', 'S');
#Uknown changes to S for Soft clipping

update probe_feature set cigar_line=replace(cigar_line, 'M', 'V');
#M will change to V (for seq&align match)

update probe_feature set cigar_line=replace(cigar_line, 'm', 'X');
#m will change to X (for seq mismatch a& align match)


INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_55_56_b.sql|probe_feature.cigar_line_format');

