
--- MAKE SURE YOU EDIT efg.sql AFTER ADDING A PATCH ! ---

-- alter cell_type.description to accommodate longer descriptions --

ALTER table cell_type change column description description varchar(80);
