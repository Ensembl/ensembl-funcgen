# patch_53_54_c.sql
#
# title: Add multispecies support
#
# description:
# Modifies coord_system and seq_region to ensure that the primary keys and structure allows for multispecies 

ALTER TABLE meta ADD KEY `meta_species_index` (`species_id`);

-- Auto-inc fields are temporarly indexed & then removed to avoid an ERROR 1075 (42000).

ALTER TABLE coord_system ADD column `species_id` int(10) DEFAULT 1;
ALTER TABLE coord_system ADD index `coord_species_idx` (`species_id`);
ALTER TABLE coord_system ADD UNIQUE index `cs_uniq_idx` (`core_coord_system_id`, `schema_build`, `species_id`);

ALTER TABLE coord_system ADD index `tmp_autoinc_idx` (`coord_system_id`);
ALTER TABLE coord_system DROP primary key;
ALTER TABLE coord_system ADD primary key (`coord_system_id`);
ALTER TABLE coord_system DROP index tmp_autoinc_idx;


ALTER TABLE seq_region ADD index `tmp_autoinc_idx` (`seq_region_id`);
ALTER TABLE seq_region DROP PRIMARY KEY;
ALTER TABLE seq_region ADD PRIMARY KEY  (`seq_region_id`);
ALTER TABLE seq_region ADD UNIQUE INDEX `sr_uniq_idx` (`name`, `schema_build`, `coord_system_id`);
ALTER TABLE seq_region drop index `tmp_autoinc_idx`;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_53_54_c.sql|multispecies');
