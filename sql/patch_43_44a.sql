
ALTER TABLE coord_system ADD `core_coord_system_id` int(10) NOT NULL;
drop index rank on coord_system;
drop index name on coord_system;
ALTER TABLE coord_system  ADD KEY `name_version_idx` (`name`, `version`);
ALTER TABLE coord_system change schema_build schema_build varchar(8) default NULL;

select "You now need to manually update the core_coord_system_ids and resolve any name version conflicts before continuing";
exit;


ALTER TABLE coord_system DROP PRIMARY KEY, ADD PRIMARY KEY(`coord_system_id`, `core_coord_system_id`, `schema_build`);

select "All name version pairs should correspond to the same(nr) coord_system_id, now update your feature and meta tables if required";
exit;

-- array.type and new key
alter table array add UNIQUE KEY   (`vendor`, `name`)
alter table array add `type` varchar(20) default NULL;
select "You need to manually update your array.type as OLIGO or PCR";



-- probe array_chip_idx
alter table probe add KEY `array_chip_idx` (`array_chip_id`);


-- probe_feature.cigar_line
-- alter table probe_feature add `cigar_line` text;
