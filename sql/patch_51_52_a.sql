# patch_51_52_a.sql
#
# title: Widen some text columns
#
# description:
# Change analysis.parameters to text, xref.description to text, coord_system.version to varchar(255)

ALTER TABLE analysis MODIFY parameters TEXT;

ALTER TABLE xref MODIFY description TEXT;

ALTER TABLE coord_system MODIFY version VARCHAR(255) DEFAULT NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_51_52_a.sql|widen_columns');


