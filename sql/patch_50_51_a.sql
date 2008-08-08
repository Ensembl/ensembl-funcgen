
--- MAKE SURE YOU EDIT efg.sql AFTER ADDING A PATCH ! ---

-- alter cell_type.description to accommodate longer descriptions --

ALTER table cell_type change column description description varchar(80);

# patch identifier
INSERT INTO meta (meta_key, meta_value) VALUES ('patch', 'patch_50_51_a.sql|cell_type_description');

