# patch_52_53_c.sql
#
# title: external_db.type ENSEMBL
#
# description:
# Add ENSEMBL to external_db.type enum

ALTER table external_db modify `type` enum('ARRAY','ALT_TRANS','MISC','LIT','PRIMARY_DB_SYNONYM', 'ENSEMBL') default NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_e.sql|external_db.type ENSEMBL');


