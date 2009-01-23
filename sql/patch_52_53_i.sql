# patch_52_53_i.sql
#
# title: probe PRIMARY KEY
#
# description:
# Alter primary key to enable probes of same name but different array_chip i.e. Illumina

ALTER table probe add UNIQUE KEY `tmp_primary` (`probe_id`,`name`, `array_chip_id`);


ALTER table probe drop primary key;

ALTER table probe add PRIMARY KEY (`probe_id`,`name`, `array_chip_id`);
ALTER table probe drop key `tmp_primary`;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_i.sql|probe PRIMARY KEY');


