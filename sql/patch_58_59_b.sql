# patch_58_59_b.sql
#
# title: add probe.description
#
# description:
# Add description field to probe able


ALTER table probe ADD `description` varchar(255) DEFAULT NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_58_59_b.sql|probe.description');


