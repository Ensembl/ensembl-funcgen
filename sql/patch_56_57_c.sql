# patch_56_57_c.sql
#
# title: unmapped_object type enum 
#
# description:
# Add array_mapping to uo.type enum

alter table unmapped_object modify `type` enum('xref','probe2transcript', 'array_mapping') NOT NULL;

#Optionally do data patch
#update unmapped_object set type='array_mapping' where type='' and (identifier='genomic' or identifier='transcript');

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_c.sql|uo.object_type_enum');


