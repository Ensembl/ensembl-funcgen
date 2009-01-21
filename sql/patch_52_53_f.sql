# patch_52_53_f.sql
#
# title: array.class 
#
# description:
# Add class field to array table to define class of array e.g. AFFY_UTR, AFFY_ST etc (This sould be used to direct parsing). Format for expression arrays would be TRANSCRIPT.

ALTER table array add column `class` varchar(20) default NULL;
ALTER table array add UNIQUE KEY `class_name_idx` (`class`, `name`);

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_f.sql|array.class');


