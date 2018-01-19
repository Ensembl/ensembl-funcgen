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

# patch_52_53_a.sql
#
# title: Update analysis_desription definition and reg feat entry?
#
# description:
# Make key UNIQUE and change displayable to boolean


ALTER table analysis_description MODIFY `displayable` BOOLEAN NOT NULL default '1';


-- Select into tmp table to remove duplciates
DROP TABLE IF EXISTS `tmp_analysis_description`;
create table `tmp_analysis_description` (
  `analysis_id` int(10) unsigned NOT NULL,
  `description` text,
  `display_label` varchar(255) default NULL,
  `displayable` BOOLEAN NOT NULL default '1',
  `web_data` text,	
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


insert ignore into tmp_analysis_description select * from analysis_description;

delete from analysis_description;
insert into analysis_description select * from tmp_analysis_description;
drop table tmp_analysis_description;

-- Finally later description for regulatory feature? We need url embedded to link directly to reg build help.

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_52_53_b.sql|redefine analysis_description');


