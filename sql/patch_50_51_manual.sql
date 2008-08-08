-- data patch for mouse
-- How was this set, we've only ever displayed H3K4me3?

--delete s from status s, status_name sn, result_set rs where rs.name='Vienna MEFf H3K4me2' and rs.result_set_id =s.table_id and s.table_name='result_set' and s.status_name_id=sn.status_name_id and sn.name='DISPLAYABLE';
--delete s from status s, status_name sn, data_set ds where ds.name='Vienna MEFf H3K4me2' and ds.data_set_id =s.table_id and s.table_name='data_set' and s.status_name_id=sn.status_name_id and sn.name='DISPLAYABLE';


-- data patch for human


--update feature_set fs, feature_type ft set fs.feature_type_id=ft.feature_type_id where ft.name='miRanda' and fs.name='miRanda miRNA';

-- end data patch
