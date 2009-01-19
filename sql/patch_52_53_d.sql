# patch_52_53_c.sql
#
# title: identity_xref ensembl xref id start end 
#
# description:
# Modify identity_xref changing 

ALTER table identity_xref change `query_identity` `xref_identity` int(5) default NULL;
ALTER table identity_xref change `target_identity` `ensembl_identity` int(5) default NULL;
ALTER TABLE identity_xref CHANGE COLUMN hit_start xref_start INT;
ALTER TABLE identity_xref CHANGE COLUMN hit_end xref_end INT;
ALTER TABLE identity_xref CHANGE COLUMN translation_start ensembl_start INT;
ALTER TABLE identity_xref CHANGE COLUMN translation_end ensembl_end INT;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_d.sql|identity_xref ensembl xref id start end');


