# patch_54_55_e.sql
#
# title: feature_type.class update
#
# description:
# Add some more feature_type.class enums remove some old and update some entries.



alter table feature_type modify `class` enum('Insulator','DNA','Regulatory Feature','Histone','RNA','Polymerase','Transcription Factor','Transcription Factor Complex','Overlap','Regulatory Motif','Region','Enhancer','Expression','Pseudo', 'Open Chromatin', 'Search Region', 'Association Locus') DEFAULT NULL;

update feature_type set class='Open Chromatin' where name='DNase1';

#What about Region/Pseudo?
#mysql> select * from feature_type where class='Region' limit 10;
#+-----------------+-------------------------+--------+----------------------------------------------------+
#| feature_type_id | name                    | class  | description                                        |
#+-----------------+-------------------------+--------+----------------------------------------------------+
#|              45 | VISTA Target            | Region | VISTA target region                                | 
#|              46 | VISTA Target - Negative | Region | Enhancer negative region identified by VISTA assay | 
#|              48 | cisRED Search Region    | Region | cisRED search region                               | 
#|          178402 | eQTL                    | Region | Expression Quantitative Trait Loci                 | 
#+-----------------+-------------------------+--------+----------------------------------------------------+


update feature_type set class='Search Region' where name in ('VISTA Target', 'VISTA Target - Negative', 'cisRED Search Region');

update feature_type set class='Association Locus' where name='eQTL';

#Now remove region/overlap
alter table feature_type modify `class` enum('Insulator','DNA','Regulatory Feature','Histone','RNA','Polymerase','Transcription Factor','Transcription Factor Complex','Regulatory Motif','Enhancer','Expression','Pseudo', 'Open Chromatin', 'Search Region', 'Association Locus') DEFAULT NULL;


INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_e.sql|feature_type.class_update');

