# patch_56_57_i.sql
#
# title: cell_type.gender
#
# description:
# Add gender colum to cell_type tables


alter table cell_type add `gender` enum('male', 'female') default NULL;


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_i.sql|cell_type.gender');


 
