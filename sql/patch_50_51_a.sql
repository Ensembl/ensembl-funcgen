
--- MAKE SURE YOU EDIT efg.sql AFTER ADDING A PATCH ! ---

-- alter cell_type.description to accommodate longer descriptions --

ALTER table cell_type change column description description varchar(80);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_50_51_a.sql|cell_type_description');


-- data patch for mouse
-- How was this set, we've only eve displayed H3K4me3?

--delete s from status s, status_name sn, result_set rs where rs.name='Vienna MEFf H3K4me2' and rs.result_set_id =s.table_id and s.table_name='result_set' and s.status_name_id=sn.status_name_id and sn.name='DISPLAYABLE';
--delete s from status s, status_name sn, data_set ds where ds.name='Vienna MEFf H3K4me2' and ds.data_set_id =s.table_id and s.table_name='data_set' and s.status_name_id=sn.status_name_id and sn.name='DISPLAYABLE';

-- end data patch



