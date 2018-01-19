-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

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


 
