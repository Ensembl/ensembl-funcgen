-- THese are patches specific to the ensembl production DB and should be done manually to ensure validity


-- data patch for mouse
-- How was this set, we've only ever displayed H3K4me3?

--delete s from status s, status_name sn, result_set rs where rs.name='Vienna MEFf H3K4me2' and rs.result_set_id =s.table_id and s.table_name='result_set' and s.status_name_id=sn.status_name_id and sn.name='DISPLAYABLE';
--delete s from status s, status_name sn, data_set ds where ds.name='Vienna MEFf H3K4me2' and ds.data_set_id =s.table_id and s.table_name='data_set' and s.status_name_id=sn.status_name_id and sn.name='DISPLAYABLE';

--delete s from status s, status_name sn, feature_set fs where fs.name='Vienna MEFf H3K4me2' and fs.feature_set_id =s.table_id and s.table_name='feature_set' and s.status_name_id=sn.status_name_id and sn.name='DISPLAYABLE';



--remove old ftype
-- delete from feature_type where name='cisRED';
--correct cisRED fset ftype
-- update feature_set fs, feature_type ft set fs.feature_type_id=ft.feature_type_id where ft.name='cisRED Motif' and fs.name='cisRED motifs';

-- data patch for human


--update feature_set fs, feature_type ft set fs.feature_type_id=ft.feature_type_id where ft.name='miRanda' and fs.name='miRanda miRNA';

--correct misRanda fset ftype

-- human and mouse  

-- update feature_set set name ='cisRED motifs' where name='cisRED group motifs'; 
-- update feature_type set description='cisRED motif' where name='cisRED Motif'; 

-- end data patch


-- Also note all cisRED data was deleted and reloaded using new group corrected definitions



