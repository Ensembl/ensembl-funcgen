# patch_52_53_c.sql
#
# title: Manual species dependent patch 
#
# description:
# Add class field to array table to define class of array e.g. AFFY_UTR, AFFY_ST etc (This sould be used to direct parsing). Format for expression arrays would be TRANSCRIPT.


SELECT "This patch is to be applied manually as it is species specific";

exit;

UPDATE external_db set db_name=replace(db_name, 'ensembl', 'homo_sapiens') where db_name like 'ensembl%';


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_manual.sql');


