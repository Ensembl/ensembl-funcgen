-- This one is for general patche
-- Can be applied multiple times with any problems
-- Put dependent patches in later patches

--Some v49 patches which were overwritten 
alter table data_set modify name varchar(100) default NULL;
alter table feature_set modify name varchar(100) default NULL;

-- Alter feature_type and add PSEUDO feature types for regulatory string
alter table feature_type modify `class` enum('Insulator','DNA','Regulatory Feature','Histone','RNA','Polymerase','Transcription Factor','Transcription Factor Complex','Overlap','Regulatory Motif','Region','Enhancer','Expression', 'Pseudo') default NULL;

