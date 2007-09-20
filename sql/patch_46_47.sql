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




-- finally remove regulatory feature from annotated_feature table
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
alter table feature_type change class class enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Overlap', 'Regulatory Motif', 'Region') default NULL; 

















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
