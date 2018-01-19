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

select "This patch should be worked through manually as there are data specific and species specific patches";
exit;




---------Need to update these in efg.sql after release!


-- meta
update meta set meta_value=48 where meta_key='schema_version';


--Ensembl only data patch human
--remove duplicate probe_features from stunnenburg set?
--delete from probe_feature where analysis_id=5;
--delete from analysis where analysis_id=5;
--remove spurious reg feat records
--delete from regulatory_feature where regulatory_feature_id=779946;
--patch ctcf result_set as displayable
--insert into status values(18,'result_set' ,2);




--remove spurious reg feat records
delete from regulatory_feature where regulatory_feature_id=779946;


--change object_xref ensembl_object_type
alter table object_xref change ensembl_object_type  ensembl_object_type ENUM('RegulatoryFeature', 'ExternalFeature') not NULL;



-- Ensembl only data patch
-- cisRED, miRanda and Vista sets, feature_types and features were deleted
-- regulatory feature type feature_type_ids were shifted, and updated in relevant tables to reset the feature_type_id to a sensible figure
--114 ->37
--115 ->38
--116 ->39
--117 ->40
--398683->41
--398684->42
--398685->43
--398686->44

--update feature_type set feature_type_id=44 where feature_type_id=398686;
--update regulatory_feature set feature_type_id=44 where feature_type_id=398686;

-- copy to tmp table to reset the autoinc counter

--create table tmp_feature_type( 
--`feature_type_id` int(10) unsigned NOT NULL auto_increment,
--   `name` varchar(40) NOT NULL,
--   `class` enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Overlap', 'Regulatory Motif', 'Region', 'Enhancer') default NULL,
--   `description`  varchar(255) default NULL,
--   PRIMARY KEY  (`feature_type_id`),
--   UNIQUE KEY `name_class_idx` (`name`, `class`)
--) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--insert into tmp_feature_type select * from feature_type;

--DROP TABLE IF EXISTS `feature_type`;
--CREATE TABLE `feature_type` (
--   `feature_type_id` int(10) unsigned NOT NULL auto_increment,
--   `name` varchar(40) NOT NULL,
--   `class` enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Overlap', 'Regulatory Motif', 'Region', 'Enhancer') default NULL,
--   `description`  varchar(255) default NULL,
--   PRIMARY KEY  (`feature_type_id`),
--   UNIQUE KEY `name_class_idx` (`name`, `class`)
--) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--insert into feature_type select * from tmp_feature_type;
--drop table tmp_feature_type;

-- added core external_dbs

-- +----------------+-----------------+------------+-----------+------------------------+------------------------+----------+-------------------------+------+-------------------+--------------------+
--| external_db_id | db_name         | db_release | status    | dbprimary_acc_linkable | display_label_linkable | priority | db_display_name         | type | secondary_db_name | secondary_db_table |
--+----------------+-----------------+------------+-----------+------------------------+------------------------+----------+-------------------------+------+-------------------+--------------------+
--|              1 | core_gene       | NULL       | KNOWNXREF |                      1 |                      0 |        5 | ensembl_core_gene       | MISC | NULL              | NULL               |
--|              2 | core_transcript | NULL       | KNOWNXREF |                      1 |                      0 |        5 | ensembl_core_transcript | MISC | NULL              | NULL               |
--+----------------+-----------------+------------+-----------+------------------------+------------------------+----------+-------------------------+------+-------------------+--------------------+


-- cisRED, miRanda and vista re-imported using new parsers




--Need to patch supporting_set, migrate type from data_set to allow multiple types of supporting set

--change channel to enum CONTROL/EXPERIMETNAL can we not NULL this also? this will default ot 1st no?


--compare mysql median to perl median for ResultFeature query

--need to vhange average row length on this and on regulatory feature!!



-- tidy up overlap feature_sets and create data_sets for them?
-- reduce size of name field in feature_set


-- should change max rows to 17000000 for reg feats to reflect max stable ids
-- also consider average row?
-- change small table primary key ids to medium int?



--add key on ec cell_type_id?
-- add enum on channel type TOTAL, EXPERIMENTAL & DUMMY? channels

-- add core regulatory tables as other_feature tables
-- regulatory_factor_coding is empty and unused? migrate to xrefs


--enum feature_type class default NULL
