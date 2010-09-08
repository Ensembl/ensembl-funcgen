# patch_59_60_d.sql
#
# title: probe_feature.cigar_line varchar
#
# description:
# Change the cigar_line field to a varchar

alter table probe_feature modify cigar_line varchar(50) default NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_59_60_d.sql|probe_feature.cigar_line_varchar');


