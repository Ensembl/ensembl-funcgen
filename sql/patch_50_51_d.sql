
--- Title: regulatory_feature.display_label

--- Description: Increase length of display_label field to accomodate longer regulatory strings

ALTER table regulatory_feature change column display_label display_label varchar(80);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_50_51_d.sql|regulatory_feature.display_label');
