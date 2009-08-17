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

