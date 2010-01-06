# patch_56_57_h.sql
#
# title: chip_seq_result_sets
#
# description:
# Alter result_set and chip channel to allow the creation of chip seq alignment
# result sets from input_sets 
#



#change table_name to enum?
#Changing table name here will prevent API working with older version of the DB
#May aswell do the input_(sub)set change at the same time.

CREATE TABLE `result_set_input` (
  `result_set_input_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `result_set_id` int(10) unsigned NOT NULL,
  `table_id` int(10) unsigned NOT NULL,
  `table_name` enum('experimental_chip', 'channel', 'input_set') default NULL,
  PRIMARY KEY (`result_set_input_id`,`result_set_id`),
  UNIQUE KEY `rset_table_idname_idx` (`result_set_id`,`table_id`,`table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


insert into result_set_input select * from chip_channel;

alter table result change `chip_channel_id` `result_set_input_id` int(10) unsigned NOT NULL;
alter table result drop key chip_channel_idx;
alter table result add key `result_set_input_idx` (`result_set_input_id`);

# Finally drop the old table
DROP table chip_channel;




# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_56_57_h.sql|chip_seq_result_set');


 
