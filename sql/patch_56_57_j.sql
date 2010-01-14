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


 
