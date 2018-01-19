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

-- create data_set and feature_set tables and entries

CREATE TABLE `feature_set` (
   `feature_set_id` int(10) unsigned NOT NULL auto_increment,
   `feature_type_id` int(10) unsigned NOT NULL,
   `analysis_id`  int(10) unsigned default NULL,
   `cell_type_id` int(10) unsigned default NULL,
   PRIMARY KEY  (`feature_set_id`),
   KEY `feature_type_idx` (`feature_type_id`)	
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `data_set` (
   `data_set_id` int(10) unsigned NOT NULL auto_increment,
   `result_set_id` int(10) unsigned default NULL,
   `feature_set_id` int(10) unsigned default NULL,
   PRIMARY KEY  (`data_set_id`, `result_set_id`, `feature_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


alter table predicted_feature add `feature_set_id` int(10) unsigned NOT NULL default '0';

-- correct feature type bodge on experimental_chip
--exp 6 & 7 ftyp_id 1 > 6
--exp 2 ftype_id 6>2




-- create feature_sets
 insert into feature_set (feature_set_id, feature_type_id, cell_type_id, analysis_id)   select distinct(ep.experiment_id), ec.feature_type_id, ec.cell_type_id, 7 from experiment_prediction ep, experimental_chip ec where ep.experiment_id=ec.experiment_id;

-- update predicted_feature and remove old columns
update predicted_feature pf, experiment_prediction ep set pf.feature_set_id=ep.experiment_id where ep.predicted_feature_id=pf.predicted_feature_id;

alter table predicted_feature drop INDEX analysis_idx;
alter table predicted_feature drop analysis_id;
alter table predicted_feature drop INDEX hit_idx;
alter table predicted_feature drop INDEX type_idx;
alter table predicted_feature drop feature_type_id;
alter table predicted_feature add index feature_set_idx(feature_set_id);

-- create data_sets for sanger

--insert into data_set(feature_set_id, result_set_id) select distinct(fs.feature_set_id), rs.result_set_id from feature_set fs, result_set rs, chip_channel cc, experimental_chip ec where ec.replicate='UNKNOWN' and ec.experimental_chip_id=cc.table_id and cc.table_name='experimental_chip' and cc.result_set_id=rs.result_set_id and fs.cell_type_id=ec.cell_type_id and fs.feature_type_id=ec.feature_type_id; 

--manually create data set for nimblegen as we are only showing one of the chipsets
--insert into data_set(feature_set_id, result_set_id) values(1,12);--human
--insert into data_set(feature_set_id, result_set_id) values(1,12);--mouse

-- make all relevant data/feature/result_sets displayable;



insert into status select experimental_chip_id, 'experimental_chip', 'DISPLAYABLE' from experimental_chip; 
--remove specific nimblegen chips.
 delete s from status s, experimental_chip ec where ec.replicate='2' and ec.experimental_chip_id=s.table_id and s.table_name ='experimental_chip' and s.state='DISPLAYABLE';


insert into status select feature_set_id, 'feature_set', 'DISPLAYABLE' from feature_set; 
insert into status select result_set_id, 'result_set', 'DISPLAYABLE' from result_set;
-- remove specific nimblegen result set
-- delete from status where table_id=13 and table_name='result_set';
