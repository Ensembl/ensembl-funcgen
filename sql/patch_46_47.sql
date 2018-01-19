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


-- meta
update meta set meta_value=47 where meta_key='schema_version';

-- add stuff here to tidy up chip_seq probe entries

-- drop coord_system_id from feature tables, now linked through seq_region tale

alter table probe_feature drop column coord_system_id;
alter table annotated_feature drop column coord_system_id;



-- add new regulatory_feature table

DROP TABLE IF EXISTS `regulatory_feature`;
CREATE TABLE `regulatory_feature` (
  `regulatory_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(1) NOT NULL default '0',	
  `display_label` varchar(60) default NULL,
  `feature_type_id`	int(10) unsigned default NULL,
  `feature_set_id`	int(10) unsigned default NULL,
  `stable_id` varchar(15) default NULL,
  PRIMARY KEY  (`regulatory_feature_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


-- Do we want a build version? Default would be schema_version
-- build may not change between version, so would have to patch table, which would indicate a build change
-- do patch to avoid having text schema_build column
-- enter build.name in meta to reflect true build? sould this be renamed schema_version?
-- now to be done in extended feature_set table

-- feature_type_id - feature type class would be regaultory feature, then have all the different types
-- display label would be dynamically set to display_label else feature type name


-- we need an nr link to cell_type
-- regulatory_feature_set or can we use feature_set?
-- this would also provide link to all contributing feature_sets
-- could remove meta record?
-- analysis_id
-- to reuse feature_set we would have to either add a feature class/type/table column 
-- then generalise feature access by dynamically linking to the table in question
-- this would mean the data set would be linking feature sets to feature sets, but of a different type.
-- this would mean we would have to have AnnotatedFeature methods and RegulatoryFeature methods
-- or can we generalise?


-- we need regulatory attribute/annotation which would map back to individual annotated_features
-- how would we map to supporting_features? (core regulatory_feature)
-- this would be used to build the feature structure dependent on it's type and whether it was focus/anchor feature
-- anchor features need recording somewhere too! data_set...no? meta? Not in attribs table due to mass redundancy
-- maybe meta, but there should be a better place for this
-- should regulatory_feature be nr for rf_id and then have classes of section of RF? Or should we split this into an exon like table, or just compute on the fly? should we have levels of display where we only have a detailed glyph once we get to a certain level of zoom? default would be simply glyph i.e. do not retrieve attrib loci and process.








-- alter feature_set table to add mutli feature type functionality
alter table feature_set add column `type` enum("annotated", "regulatory", "supporting") default NULL;
update feature_set set type='annotated';

-- ensembl specific data patch
-- migrate v45 regulatory features
--insert into regulatory_feature select NULL, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, display_label, NULL, feature_set_id from annotated_feature where feature_set_id=27;

--
--update regulatory_feature set stable_id=replace(display_label, ":%", "")
--didn't work, cropped by changing field width

-- crop stable id from binary string
update regulatory_feature set display_label=substring(display_label from 17);

-- ensembl specific data patch
-- point reg feature_set at new table, but do not delete old annotated features just yet
-- update feature_set set type='regulatory' where feature_set_id=27; 


-- update regulatory feature_types
insert into feature_type values(NULL, 'Gene Associated', 'REGULATORY FEATURE', 'Gene like regulatory feature'); 
insert into feature_type values(NULL, 'Promoter Associated', 'REGULATORY FEATURE', 'Promoter like regulatory feature'); 
insert into feature_type values(NULL, 'Non-Gene Associated', 'REGULATORY FEATURE', 'Non-Gene like regulatory feature'); 
insert into feature_type values(NULL, 'Unclassified', 'REGULATORY FEATURE', 'Unclassified regulatory feature'); 


-- Table structure for table `supporting_set`

DROP TABLE IF EXISTS `supporting_set`;
CREATE TABLE `supporting_set` (
   `data_set_id` int(10) unsigned NOT NULL default '0',
   `supporting_set_id` int(10) unsigned NOT NULL default '0',
   PRIMARY KEY  (`data_set_id`, `supporting_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- primary key will always be unique as we will never 
-- have mixed supporting_set type e.g.(result/feature) in the same data_set
-- hence no possibilty of getting same id from different tables in same data_set


alter table data_set add column `supporting_set_type` enum("result", "feature") default NULL; 
-- can't have not NULL here or will default to result!  Will this default to result on other than enum val?

-- populate data_set supporting_set_type
update data_set set supporting_set_type='result';

-- ensembl specific data patches 
--update data_set set supporting_set_type='feature' where name like "Regulatory_%";

-- populate the supporting_set table
insert into supporting_set select data_set_id, result_set_id from data_set;
--tidy data_set_memebr table where we only have feature_sets
delete from supporting_set where supporting_set_id=0;


-- remove old result_set_id column from data_set
alter table data_set drop column result_set_id; 


-- alter keys on data_set
alter table data_set add KEY `supporting_type_idx` (`supporting_set_type`);



-- now data_set can handled feature to feature sets
-- add the regulatory_attribute table to capture the specific feature overlaps/attributes

-- add new regulatory_attribute table

DROP TABLE IF EXISTS `regulatory_attribute`;
CREATE TABLE `regulatory_attribute` (
  `regulatory_feature_id` int(10) unsigned NOT NULL default '0',
  `attribute_feature_id` int(10) unsigned NOT NULL default '0',
  `attribute_feature_table` enum('annotated', 'supporting') default NULL,
  PRIMARY KEY  (`regulatory_feature_id`, `attribute_feature_table`, `attribute_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=17;



-- do we need key on feature_type and or feature_id?



-- check we have migrated the regulatory_feature before this
-- finally remove regulatory feature from annotated_feature table???
delete from annotated_feature where feature_set_id = 27;


-- patched regulatory_feature feature_typ and stable_id using patch_reg_feature.pl
-- no alter column defs for stable_id
alter table regulatory_feature change column stable_id `stable_id` mediumint(8) unsigned default NULL;
-- allows for ~16 millions regulatory features
--add index
alter table regulatory_feature add  KEY `stable_id_idx` (`stable_id`);

-- v47 build specific data tidy up
-- delete from annotated_feature where feature_set_id >56;
-- update feature_set set type='regulatory' where feature_set_id >56;



-- tweak result table
-- No need for avg length or result_id as primary key
-- need probe_id as first in primary key
-- need chip_channel_idx for chip level methods, i.e. norm

CREATE TABLE `tmp_result` (
   `result_id` int(10) unsigned NOT NULL auto_increment,
   `probe_id` int(10) unsigned default NULL,
   `score` double default NULL,
   `chip_channel_id` int(10) unsigned NOT NULL,
   `X` int(4) unsigned default NULL,
   `Y` int(4) unsigned default NULL,
   PRIMARY KEY  (`result_id`),
   KEY `probe_idx` (`probe_id`),
   KEY `chip_channel_idx` (`chip_channel_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;


--do we really need chip_channel_idx, or will this just be a romp through the table anyway.
--probably can get rid when we federate, will reduce the size of the DB.

insert into tmp_result select * from result;


DROP TABLE IF EXISTS `result`;
CREATE TABLE `result` (
   `result_id` int(10) unsigned NOT NULL auto_increment,
   `probe_id` int(10) unsigned default NULL,
   `score` double default NULL,
   `chip_channel_id` int(10) unsigned NOT NULL,
   `X` smallint(4) unsigned default NULL,
   `Y` smallint(4) unsigned default NULL,
   PRIMARY KEY  (`result_id`),
   KEY `probe_idx` (`probe_id`),
   KEY `chip_channel_idx` (`chip_channel_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

insert into result select * from tmp_result;
DROP TABLE IF EXISTS `tmp_result`;



--
-- Table structure for table `experimental_set`
-- experimental_file?

-- This is to accomodate any direct import of features from a pre-processed experiment
-- e.g. a short reads experiment
-- It by passes all the array_chip/experimental_chip/channel/result level information and ties directly into a feature/data_set
-- we need to be able to record the status of each set member/file individually for recovery purposes?
-- No because we can't roll bak on a file basis anyway due to lack of chip_channel like ids in the annotated_feature table
-- Would this be restricted to one feature per experiment?
-- If not need to have a experimental_set_id
-- best to do this anyway so we are in line with chip experiemtns in allowing multiple feature/cell_type per experimenht
-- and then we also have an experimental_set_idin the data_set table

-- should we intercede 'file' into ther table names to make it more clear? 

DROP TABLE IF EXISTS `experimental_set`;
CREATE TABLE `experimental_set` (
   `experimental_set_id` int(10) unsigned NOT NULL auto_increment,
   `experiment_id` int(10) unsigned default NULL, --
   `feature_type_id` int(10) unsigned default NULL,
   `cell_type_id` int(10) unsigned default NULL,
   `format` varchar(20) default NULL,
   `vendor` varchar(40) default NULL,
   `name` varchar(40) not NULL default '0',
   PRIMARY KEY  (`experimental_set_id`),
   UNIQUE KEY `name_idx` (`name`),
   KEY `experiment_idx` (`experiment_id`),
   KEY `feature_type_idx` (`feature_type_id`),
   KEY `cell_type_idx` (`cell_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=30;


-- do we want a name field? 40 to coomodate exp_name plus BR /TR notation etc.
-- do we want a type field? 'SHORT_READ' or ????
-- format is type of short read, platform name? 
-- now where do we put the experimental type? Format? Vendor?
-- do we need an auxilliary table here akin to supporting set to remove redundancy of cell_type, feature_type
-- Do we need replicates if we are peak calling outside the DB?
-- keys on vendor/format?

--
-- Table structure for table `experimental_subset`
--

-- represents a file from a subset
-- mainly used for tracking import and recovery 

DROP TABLE IF EXISTS `experimental_subset`;
CREATE TABLE `experimental_subset` (
   `experimental_subset_id` int(10) unsigned NOT NULL auto_increment,
   `experimental_set_id` int(10) unsigned NOT NULL default '0',
   `name` varchar(30) NOT NULL default '0', -- filename?	
   PRIMARY KEY  (`experimental_subset_id`), 
   UNIQUE KEY `set_name_dx` (`experimental_set_id`, `name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=30;


-- remove experimental_set_id from key to force uniqueness of file?
-- No file names may be same across expeirments, this allows addition of same sub set to different sets
-- whilst making it unique whtin a set
-- no we can't have multi experiment sets at present, still migh have same filename tho


--uppdate data_set
alter table data_set change `supporting_set_type` `supporting_set_type` enum('result', 'feature', 'experimental') default NULL;

--Ensembl data pach only
--tidy up old overlap features/set/types
delete from annotated_feature where feature_set_id in(56, 55);
delete from feature_type where name='CTCF:CTCF:DNase1:H2AZ:H2BK5me1:H3K27me1:';
delete from feature_set where feature_set_id in(56, 55);
update feature_type set name='Wiggle_H3K4me3_focus' where feature_type_id=17;
update feature_type set name='Nessie_NG_STD_2_ctcf_ren_BR1_focus' where feature_type_id=16;
update feature_type set name='GM06990_DNASE_IMPORT_focus' where feature_type_id=18;


-- Tidy up sql discrepencies

alter table array change array_id `array_id` int(10) unsigned NOT NULL auto_increment;
alter table array_chip change array_chip_id `array_chip_id` int(10) unsigned NOT NULL auto_increment;
alter table array_chip change array_id `array_id` int(10) unsigned NOT NULL;
alter table array_chip drop key design_idx;
alter table array_chip drop key array_idx;
alter table array_chip add unique key `array_design_idx` (`array_id`, `design_id`);
alter table channel change channel_id `channel_id` int(10) unsigned NOT NULL auto_increment;
alter table channel change experimental_chip_id `experimental_chip_id` int(11) unsigned default NULL;

alter table meta_coord change coord_system_id `coord_system_id` int(10) NOT NULL;
alter table coord_system change coord_system_id `coord_system_id` int(10) NOT NULL auto_increment;

alter table experimental_chip change experimental_chip_id `experimental_chip_id` int(10) unsigned NOT NULL auto_increment;
alter table experimental_chip change experiment_id `experiment_id` int(10) unsigned NOT NULL;
alter table experimental_chip change array_chip_id `array_chip_id` int(10) unsigned NOT NULL;
alter table experimental_chip drop key chip_idx;
alter table experimental_chip add KEY `feature_type_idx` (`feature_type_id`);
alter table experimental_chip add KEY `unique_id_idx` (`unique_id`);

alter table experiment change experiment_id `experiment_id` int(10) unsigned NOT NULL auto_increment;


alter table feature_type change feature_type_id `feature_type_id` int(10) unsigned NOT NULL auto_increment;
alter table feature_type drop key feature_type_name_idx;
alter table feature_type add UNIQUE KEY `name_class_idx` (`name`, `class`);


alter table meta change meta_id `meta_id` int(10) NOT NULL auto_increment;
alter table status change table_id `table_id` int(10) unsigned NOT NULL default '0';
alter table experimental_variable change table_id `table_id` int(10) unsigned NOT NULL default '0';

alter table regulatory_attribute change `attribute_feature_table` `attribute_feature_table` enum('annotated', 'external') NOT NULL default 'annotated';
alter table data_set change  supporting_set_type `supporting_set_type` enum('result', 'feature', 'experimental', 'external') default NULL;
alter table feature_set change type `type` enum("annotated", "regulatory", "external") default NULL;


--add xref schema


--xref stuff

DROP TABLE IF EXISTS `object_xref`;
CREATE TABLE object_xref (
  object_xref_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  ensembl_id                  INT(10) UNSIGNED NOT NULL, 
  ensembl_object_type         ENUM('regulatory_feature', 'external_feature')
                              not NULL,
  xref_id                     INT UNSIGNED NOT NULL,
  linkage_annotation          VARCHAR(255) DEFAULT NULL,
  UNIQUE (ensembl_object_type, ensembl_id, xref_id),
  KEY oxref_idx (object_xref_id, xref_id, ensembl_object_type, ensembl_id),
  KEY xref_idx (xref_id, ensembl_object_type)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=40;


--regulatory_feature or RegulatoryModule or/and RegulatoryRegion?
--we are going to have to do the vega trick here of loading all the regulatory features as an external_db to enable xrefs to core.
--ensembl_id could be eFG dbID or core stable_id

DROP TABLE IF EXISTS identity_xref;
CREATE TABLE identity_xref (
  object_xref_id          INT(10) UNSIGNED NOT NULL,
  query_identity 	  INT(5),
  target_identity         INT(5),
  hit_start               INT,
  hit_end                 INT,
  translation_start       INT,
  translation_end         INT,
  cigar_line              TEXT, 
  score                   DOUBLE,
  evalue                  DOUBLE,
  analysis_id             SMALLINT UNSIGNED NOT NULL,
  PRIMARY KEY (object_xref_id),
  KEY analysis_idx (analysis_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS xref;
CREATE TABLE xref (
   xref_id 		      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
   external_db_id             SMALLINT UNSIGNED NOT NULL,
   dbprimary_acc              VARCHAR(40) NOT NULL,
   display_label              VARCHAR(128) NOT NULL,
   version                    VARCHAR(10) DEFAULT '0' NOT NULL,
   description                VARCHAR(255),
   info_type                  ENUM('PROJECTION', 'MISC', 'DEPENDENT', 'DIRECT', 'SEQUENCE_MATCH', 'INFERRED_PAIR', 'PROBE', 'UNMAPPED') not NULL,
   info_text                  VARCHAR(255),
   PRIMARY KEY (xref_id),
   UNIQUE KEY id_index (dbprimary_acc, external_db_id, info_type, info_text),
   KEY display_index (display_label),
   KEY info_type_idx (info_type)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=100;



--  Table structure for table 'external_synonym'

DROP TABLE IF EXISTS external_synonym;
CREATE TABLE external_synonym (
  xref_id                     INT(10) UNSIGNED NOT NULL,
  synonym                     VARCHAR(40) NOT NULL, 
  PRIMARY KEY (xref_id, synonym),
  KEY name_index (synonym)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=20;



-- Table structure for table 'external_db' 

DROP TABLE IF EXISTS external_db;
CREATE TABLE external_db (
  external_db_id 	          SMALLINT(5) UNSIGNED NOT NULL auto_increment,
  db_name                     VARCHAR(28) NOT NULL,
  db_release                  VARCHAR(255),
  status                      ENUM('KNOWNXREF','KNOWN','XREF','PRED','ORTH', 'PSEUDO') NOT NULL,
  dbprimary_acc_linkable      BOOLEAN DEFAULT 1 NOT NULL,
  display_label_linkable      BOOLEAN DEFAULT 0 NOT NULL,
  priority                    INT NOT NULL,
  db_display_name             VARCHAR(255),
  type                        ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM'),
  secondary_db_name           VARCHAR(255) DEFAULT NULL,
  secondary_db_table          VARCHAR(255) DEFAULT NULL,
  PRIMARY KEY (external_db_id) 
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=80;



-- Table structure for table 'go_xref'

DROP TABLE if EXISTS go_xref;
CREATE TABLE go_xref (
  object_xref_id          INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  linkage_type            ENUM('IC', 'IDA', 'IEA', 'IEP', 'IGI', 'IMP', 
		               'IPI', 'ISS', 'NAS', 'ND', 'TAS', 'NR', 'RCA') NOT NULL,
  source_xref_id          INT(10) UNSIGNED DEFAULT NULL,
  KEY (object_xref_id),
  KEY (source_xref_id),
  UNIQUE (object_xref_id, source_xref_id, linkage_type)
)  ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- This is just an empty to table to avoid having to rework all the core sql and API to accomodate eFG specific xref schema


-- back to the external_feature stuff

--prepare feature_type table to receive new feature classes

update feature_type set class='Insulator' where class='INSULATOR';
update feature_type set class ='Histone' where class ='Histone';
update feature_type set class ='Regulatory Feature' where class ='REGULATORY FEATURE';
update feature_type set class ='Overlap' where class ='OVERLAP';
update feature_type set class ='Polymerase' where class ='POLYMERASE';
update feature_type set class ='Transcription Factor' where class ='TRANSCRIPTION FACTOR';

alter table feature_type change class class enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Overlap', 'Regulatory Motif', 'Region', 'Enhancer') default NULL; 



DROP TABLE IF EXISTS `external_feature`;
CREATE TABLE `external_feature` (
  `external_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(1) NOT NULL default '0',	
  `display_label` varchar(60) default NULL,
  `feature_type_id`	int(10) unsigned default NULL,
  `feature_set_id` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`external_feature_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=80;




insert into feature_type values(NULL, 'cisRED Search Region', 'Region', 'cisRED search region');
insert into feature_type values(NULL, 'cisRED', 'Regulatory Motif', 'cisRED group motif set');
insert into feature_type values(NULL, 'miRanda', 'RNA', 'miRanda microRNA set');


insert into analysis(name) values('miRanda');
insert into analysis_description select analysis_id, 'miRanda microRNA target prediction (http://cbio.mskcc.org/mirnaviewer/)', 'miRanda', 0, NULL from analysis where logic_name='miRanda';


insert into analysis(logic_name) values('cisRED');
insert into analysis_description select analysis_id, 'cisRED motif search (www.cisred.org)', 'cisRED', 0, NULL from analysis where logic_name='cisRED';



--Now the data patch for human only
--this is an ensembl only internal patch

--dump the following tables from v46 core
--regulatory_feature
--regulatory_factor
--regulatory_factor_coding
--regulatory_search_region
--edit the regulatory_feature dump to regulatory_feature_core(or copy and patch table before dump)
--to avoid overwriting the eFG regulatory_feature_table

--import all above tables

-- we have 4 distinct sets from 3 sources

--cisRED atomic patterns(feature_types) and group motifs (external_features)
--cisRED search regions, 1 type with multiple features
--miRanda microRNAs(feature_type) and alignments(external_feature)
--http://enhancer.lbl.gov/ enhancer positive regions and enhancer negative regions.

--each has a generic feature_type to define the feature_set, and (apart from the cisRED regions) at least 2 or many more feature_type to represent alignments/motif hits or +ve or -ve regions.

-- we need to merge enhancer analyses into one, so we can have one set giving different feature_types at the external_feature level.
-- also merge cisRED analyses into one and represent search regions and motif simply as sifferent typs in two different sets

--is the enhancer set just Mouse data?

-- let's deal with the miRanda first
-- first let's create set with all the factor orphans removed i.e. factor which don't have features
-- there are numerous factors which don't have mapping, most of which are present 3-10 times with different ids but the same name


CREATE TABLE `regulatory_factor_nr` (
  `regulatory_factor_id` int(10) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `type` varchar(255) default NULL,
  PRIMARY KEY  (`regulatory_factor_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into regulatory_factor_nr select rf.* from regulatory_factor rf, regulatory_feature_core rfc where rf.regulatory_factor_id=rfc.regulatory_factor_id group by rf.regulatory_factor_id;


-- now migrate the miRanda factors to feature_types
-- we we need to clean up the non-hsa ones? These are already removed by the feature only select above
insert into feature_type select regulatory_factor_id, name, 'RNA', 'miRanda miRNA' from regulatory_factor_nr where type='miRNA_target';

--create the feature_set
insert into feature_set select NULL, ft.feature_type_id, a.analysis_id, NULL, 'miRanda miRNA', 'external' from feature_type ft, analysis a where ft.name='miRanda' and a.logic_name='miRanda';

--now migrate the features
insert into external_feature select NULL, rfc.seq_region_id, rfc.seq_region_start, rfc.seq_region_end, rfc.seq_region_strand, rfc.name, ft.feature_type_id, fs.feature_set_id from regulatory_feature_core rfc, feature_type ft, feature_set fs, regulatory_factor_nr rfn where fs.name='miRanda miRNA' and ft.name=rfn.name and rfc.regulatory_factor_id=rfn.regulatory_factor_id;


-- now the cisRED search set
insert into feature_set select NULL, ft.feature_type_id, a.analysis_id, NULL, 'cisRED search regions', 'external' from feature_type ft, analysis a where ft.name='cisRED Search Region' and a.logic_name='cisRED';


--migrate the features from regulatory_search_region
--this has been redone with v47 regions which were mapped to NCBI36

insert into external_feature select NULL, rsr.seq_region_id, rsr.seq_region_start, rsr.seq_region_end, rsr.seq_region_strand, rsr.name, fs.feature_type_id, fs.feature_set_id from regulatory_search_region rsr, feature_set fs where fs.name='cisRED Search Regions';

--we will have to patch the xrefs later so don't delete this table!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


-- now the cisRED motifs
insert into feature_set select NULL, ft.feature_type_id, a.analysis_id, NULL, 'cisRED group motifs', 'external' from feature_type ft, analysis a where ft.name='cisRED' and a.logic_name='cisRED';


--were going to do the select on is null, so need to change this for the empty enhancer factor record
update regulatory_factor set name='Vista Enhancer', type='Enhancer' where name ='';
--probably need to rename this

insert into feature_type select regulatory_factor_id, name, 'Regulatory Motif', 'cisRED group motif' from regulatory_factor_nr where type is null;

insert into external_feature select NULL, rfc.seq_region_id, rfc.seq_region_start, rfc.seq_region_end, rfc.seq_region_strand, rfc.name, ft.feature_type_id, fs.feature_set_id from regulatory_feature_core rfc, feature_type ft, feature_set fs, regulatory_factor_nr rfn where fs.name='cisRED group motifs' and ft.name=rfn.name and rfc.regulatory_factor_id=rfn.regulatory_factor_id;

--now finally the enhancer set
insert into feature_type values(NULL, 'VISTA Target', 'Region', 'VISTA target region');
insert into feature_type values(NULL, 'VISTA Enhancer', 'Enhancer', 'Enhancer identified by positive VISTA assay');
-- do we need this last one as? It is implicit that there is no enhancer
-- do we need to change the name??
insert into feature_type values(NULL, 'VISTA Target - Negative', 'Region', 'Enhancer negative region identified by VISTA assay');

insert into analysis(logic_name) values('VISTA');
insert into analysis_description select analysis_id, 'VISTA Enhancer Assay (http://enhancer.lbl.gov/)', 'VISTA', 0, NULL from analysis where logic_name='VISTA';


insert into feature_set select NULL, ft.feature_type_id, a.analysis_id, NULL, 'VISTA enhancer set', 'external' from feature_type ft, analysis a where ft.name='VISTA Target' and a.logic_name='VISTA';

-- +ve features frist analysis_id 5210
insert into external_feature select NULL, rfc.seq_region_id, rfc.seq_region_start, rfc.seq_region_end, rfc.seq_region_strand, rfc.name, ft.feature_type_id, fs.feature_set_id from regulatory_feature_core rfc, feature_type ft, feature_set fs, regulatory_factor_nr rfn where ft.name='VISTA Enhancer'  and fs.name='VISTA enhancer set' and rfc.regulatory_factor_id=rfn.regulatory_factor_id and rfc.analysis_id=5210;

-- -ve features frist analysis_id 5211
insert into external_feature select NULL, rfc.seq_region_id, rfc.seq_region_start, rfc.seq_region_end, rfc.seq_region_strand, rfc.name, ft.feature_type_id, fs.feature_set_id from regulatory_feature_core rfc, feature_type ft, feature_set fs, regulatory_factor_nr rfn where ft.name='VISTA Target - Negative'  and fs.name='VISTA enhancer set' and rfc.regulatory_factor_id=rfn.regulatory_factor_id and rfc.analysis_id=5211;


-- now we can drop some tables
drop table regulatory_factor;
drop table regulatory_factor_nr;
drop table regulatory_feature_core;
drop table regulatory_factor_coding;

-- oops we had all the miRNAs in the cisRED set :|
delete ef from feature_type ft, external_feature ef where ef.feature_type_id=ft.feature_type_id and ft.class ='RNA' and ef.feature_set_id=61;
delete from feature_type where feature_type_id=375975;
-- and the enhancer set..doh!



-- darn it, forgot about efg seq_region_ids
--need to convert to core sr_ids in external_feature table to efg sr_ids
-- they are all on cs_id 1 and we no none match between the two sets as this returns nothing
-- select distinct(coord_system_id) from seq_region sr, external_feature ef where ef.seq_region_id=sr.seq_region_id;
--this means we can do a direct update without worrying about clashing sr_ids between core annd efg.
update external_feature ef, seq_region sr set ef.seq_region_id=sr.seq_region_id where ef.seq_region_id=sr.core_seq_region_id and sr.coord_system_id=1;

--and finally meta_coord
select (seq_region_end - seq_region_start + 1) as maxlength from external_feature order by maxlength desc limit 5;
insert into meta_coord values('external_feature', 1, 131013);

-- cisRED search regions are still all on NCBI35 coords, but in a v36 seq_region_table, i.e. 35 is not the default build.
-- we need to import this seq_region info.

-- we still have the search_region table as we want to populate the xrefs

-- This is presently unfinished as we need to manually store all 35 seq_regions as toplevel will not return any.





-- Ensembl only data patch
-- Migrate all ChIP-Seq imports to the experimental_set table

insert into experimental_set select NULL, ec.experiment_id, ec.feature_type_id, ec.cell_type_id, 'SEQUENCING', 'SOLEXA', ec.unique_id from experimental_chip ec where ec.experiment_id>=13;


update data_set ds, supporting_set ss, experimental_set es set ds.supporting_set_type='experimental', ss.supporting_set_id=es.experimental_set_id where ds.name=es.name and ds.data_set_id=ss.data_set_id;

--CD4_DNASE result_set and have been duplicated for some reason, no feature set for duplciation so remove this set first;
delete ds, ss from data_set ds, supporting_set ss where ds.name='CD4_DNASE' and ds.data_set_id=ss.data_set_id;

-- Now update names in data_set and experimental_set
update data_set set name=replace(name, '_IMPORT', '') where data_set_id in (20,21);
update experimental_set set name='CD4_DNASE' where name ='CD4_parzen_02';
update experimental_set set name='GM06990_DNASE' where name='GM06990_parzen_0115';

--Now redo experimental_set patch for these
update data_set ds, supporting_set ss, experimental_set es set ds.supporting_set_type='experimental', ss.supporting_set_id=es.experimental_set_id where ds.name=es.name and ds.data_set_id=ss.data_set_id and ds.data_set_id in(20,21);

--Now tidy up all array/chip/probe/probe_feature level records for these experiments

--results first 
-- create tmp table to aid delete
DROP TABLE IF EXISTS `tmp_ccids`;
CREATE TABLE `tmp_ccids` (
  `chip_channel_id` int(10) unsigned NOT NULl,
  PRIMARY KEY  (`chip_channel_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;



insert into tmp_ccids select cc.chip_channel_id from chip_channel cc, experimental_chip ec where cc.table_name='experimental_chip' and cc.table_id=ec.experimental_chip_id and ec.experiment_id>=13;

insert into tmp_ccids select cc.chip_channel_id from chip_channel cc, experimental_chip ec, channel c where cc.table_name='channel' and cc.table_id=c.channel_id and c.experimental_chip_id=ec.experimental_chip_id and ec.experiment_id>=13;

delete r from result r, tmp_ccids tc where r.chip_channel_id=tc.chip_channel_id;
--Query OK, 24205493 rows affected (8 min 45.35 sec)

DROP TABLE IF EXISTS `tmp_ccids`;

-- now delete the channel and experimental_chip records
delete c from channel c, experimental_chip ec where c.experimental_chip_id=ec.experimental_chip_id and ec.experiment_id>=13;
delete from experimental_chip where experiment_id>=13;

--now delete the pseudo probes/features
--do another tmp table to aid the delete
DROP TABLE IF EXISTS `tmp_probe_ids`;
CREATE TABLE `tmp_probe_ids` (
  `probe_id` int(10) unsigned NOT NULL,
  PRIMARY KEY  (`probe_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

insert  into tmp_probe_ids select probe_id from probe where array_chip_id>=44;
--Query OK, 24869910 rows affected (3 min 6.35 sec)
delete p from probe p, tpi_probe_ids tpi where p.probe_id=tpi.probe_id;
--Query OK, 24869910 rows affected (20 min 53.63 sec)

delete pf from probe_feature pf, tmp_probe_ids tpi where pf.probe_id=tpi.probe_id;
--Query OK, 24869908 rows affected (10 min 54.17 sec)
-- array/chips and result sets/chip_channel_ids we forgot
delete from array_chip where array_chip_id>=44;
delete from array where array_id>=3;

delete cc from chip_channel cc, result_set rs where cc.result_set_id=rs.result_set_id and rs.result_set_id>=23;
delete from result_set where result_set_id>=23;


--That's all the external_sets cleaned up, just need to add the filenames for each set as subsets.
-- or is it?  What about the wiggle sets?


--tweak RegFeat entry in feature_type table
update feature_type set name='RegulatoryFeature', description='Regulatory Feature' where name='RegulatoryFeatures';

alter table feature_type change name `name` varchar(40) NOT NULL;







-- Data patch for overlapping old reg feats, this is to enable stable id mapping code to work


--The first features in the overlapping pairs

mysql> select * from regulatory_feature rf inner join regulatory_feature rf1 on rf1.seq_region_start=rf.seq_region_end and rf.seq_region_id=rf1.seq_region_id and rf1.feature_set_id=rf.feature_set_id and rf.feature_set_id=27;
+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+
| regulatory_feature_id | seq_region_id | seq_region_start | seq_region_end | seq_region_strand | display_label | feature_type_id | feature_set_id | stable_id | regulatory_feature_id | seq_region_id | seq_region_start | seq_region_end | seq_region_strand | display_label | feature_type_id | feature_set_id | stable_id |
+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+
|                 23522 |           257 |         12141408 |       12148061 |                 0 | 0011111100    |             117 |             27 |     23522 |                 23523 |           257 |         12148061 |       12149948 |                 0 | 0011011100    |             117 |             27 |     23523 |
|                 23541 |           257 |         12331307 |       12332398 |                 0 | 0011011100    |             117 |             27 |     23541 |                 23542 |           257 |         12332398 |       12335452 |                 0 | 0011111100    |             117 |             27 |     23542 |
|                 52037 |            29 |         45003728 |       45005056 |                 0 | 0011011100    |             117 |             27 |     52037 |                 52038 |            29 |         45005056 |       45008548 |                 0 | 0011011100    |             117 |             27 |     52038 |
|                 63403 |            89 |        171237984 |      171240084 |                 0 | 0100011010    |             117 |             27 |     63403 |                 63404 |            89 |        171240084 |      171241746 |                 0 | 0000111010    |             117 |             27 |     63404 |
|                 79529 |            41 |        141476862 |      141483419 |                 0 | 0011111100    |             117 |             27 |     79529 |                 79530 |            41 |        141483419 |      141484341 |                 0 | 0011111100    |             117 |             27 |     79530 |
|                104835 |            25 |        127359713 |      127362848 |                 0 | 0011011100    |             117 |             27 |    104835 |                104836 |            25 |        127362848 |      127364726 |                 0 | 0011011100    |             117 |             27 |    104836 |
|                104983 |            25 |        128777241 |      128778840 |                 0 | 0011011100    |             117 |             27 |    104983 |                104984 |            25 |        128778840 |      128782578 |                 0 | 0011011100    |             117 |             27 |    104984 |
|                105122 |            25 |        130477548 |      130480838 |                 0 | 0011111100    |             117 |             27 |    105122 |                105123 |            25 |        130480838 |      130482511 |                 0 | 0011111100    |             117 |             27 |    105123 |
+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+
8 rows in set (0.84 sec)



-- THe second feature in the overlapping pairs

select * from regulatory_feature rf inner join regulatory_feature rf1 on rf.seq_region_start=rf1.seq_region_end and rf.seq_region_id=rf1.seq_region_id and rf1.feature_set_id=rf.feature_set_id and rf.feature_set_id=27;



+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+
| regulatory_feature_id | seq_region_id | seq_region_start | seq_region_end | seq_region_strand | display_label | feature_type_id | feature_set_id | stable_id | regulatory_feature_id | seq_region_id | seq_region_start | seq_region_end | seq_region_strand | display_label | feature_type_id | feature_set_id | stable_id |
+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+
|                 23523 |           257 |         12148061 |       12149948 |                 0 | 0011011100    |             117 |             27 |     23523 |                 23522 |           257 |         12141408 |       12148061 |                 0 | 0011111100    |             117 |             27 |     23522 |
|                 23542 |           257 |         12332398 |       12335452 |                 0 | 0011111100    |             117 |             27 |     23542 |                 23541 |           257 |         12331307 |       12332398 |                 0 | 0011011100    |             117 |             27 |     23541 |
|                 52038 |            29 |         45005056 |       45008548 |                 0 | 0011011100    |             117 |             27 |     52038 |                 52037 |            29 |         45003728 |       45005056 |                 0 | 0011011100    |             117 |             27 |     52037 |
|                 63404 |            89 |        171240084 |      171241746 |                 0 | 0000111010    |             117 |             27 |     63404 |                 63403 |            89 |        171237984 |      171240084 |                 0 | 0100011010    |             117 |             27 |     63403 |
|                 79530 |            41 |        141483419 |      141484341 |                 0 | 0011111100    |             117 |             27 |     79530 |                 79529 |            41 |        141476862 |      141483419 |                 0 | 0011111100    |             117 |             27 |     79529 |
|                104836 |            25 |        127362848 |      127364726 |                 0 | 0011011100    |             117 |             27 |    104836 |                104835 |            25 |        127359713 |      127362848 |                 0 | 0011011100    |             117 |             27 |    104835 |
|                104984 |            25 |        128778840 |      128782578 |                 0 | 0011011100    |             117 |             27 |    104984 |                104983 |            25 |        128777241 |      128778840 |                 0 | 0011011100    |             117 |             27 |    104983 |
|                105123 |            25 |        130480838 |      130482511 |                 0 | 0011111100    |             117 |             27 |    105123 |                105122 |            25 |        130477548 |      130480838 |                 0 | 0011111100    |             117 |             27 |    105122 |
+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+-----------------------+---------------+------------------+----------------+-------------------+---------------+-----------------+----------------+-----------+
8 rows in set (0.00 sec)

--Now tweak the 3'/2nd feature to have a start value of srtart+1

--update regulatory_feature set seq_region_start=(seq_region_start+1) where regulatory_feature_id=(select rf.regulatory_feature_id from regulatory_feature rf inner join regulatory_feature rf1 on rf.seq_region_start=rf1.seq_region_end and rf.seq_region_id=rf1.seq_region_id and rf1.feature_set_id=rf.feature_set_id and rf.feature_set_id=27);

--will this work with the nested select on an update?
--need to double nest, to get work around just do in for simplicity
update regulatory_feature set seq_region_start=(seq_region_start+1) where regulatory_feature_id IN (23523, 23542,52038,  63404, 79530, 104836, 104984, 105123);


-- Now patch/remove the old overlap sets in human
delete from annotated_feature where feature_set_id in(24,25,26);
delete from feature_set where feature_set_id in(24,25,26);


--General patch! i.e. non-data patch
alter table feature_set drop index name_idx;
alter table feature_set change name `name` varchar(40) default NULL;
alter table feature_set add UNIQUE KEY `name_idx` (`name`);


-- Remove search region table
-- These contain the xrefs, which we have not yet loaded, but we will just reload them all next release 
-- when we use the import script instead of a patch
-- possible bug in migration as vista element  count doesn't match that in core
DROP table if exists regulatory_search_region;



-- mouse only meta_coord data patch
update meta_coord set max_length=50 where table_name='probe_feature';

--human data patch
--remove old CS level imported status entries
delete s, sn  from status s, status_name sn where s.status_name_id=sn.status_name_id and sn.name='IMPORTED_CS_17';

--need to vhange average row length on this and on regulatory feature!!



-- tidy up overlap feature_sets and create data_sets for them
-- reduce size of name field in feature_set


-- should change max rows to 17000000 for reg feats to reflect max stable ids
-- also consider average row?
-- change small table primary key ids to medium int?



--add key on ec cell_type_id?
-- add enum on channel type TOTAL, EXPERIMENTAL & DUMMY? channels

-- add core regulatory tables as other_feature tables
-- regulatory_factor_coding is empty and unused? migrate to xrefs


--enum feature_type class default NULL

-- correct data_set.supporting_set_type for Wiggle_H* and H3ac_claes data
-- update data_set set supporting_set_type='experimental' where data_set_id >= 12 AND data_set_id <= 19;
