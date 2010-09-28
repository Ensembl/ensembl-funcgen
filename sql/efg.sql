
# Ensembl table definitions
# 
# Note that more information about each table can be found in
# ensembl/docs/schema_description/

# Conventions:
#  - use lower case and underscores
#  - internal ids are integers named tablename_id
#  - same name is given in foreign key relations

# DO NOT use '---' as a comment, this breaks in macosx mysql?!
# This also generally applies to all  --[^ ] !?


-- CORE TABLES --

--
-- Table structure for table `analysis`
--

DROP TABLE IF EXISTS `analysis`;
CREATE TABLE `analysis` (
  `analysis_id` smallint(5) unsigned NOT NULL auto_increment,
  `created` datetime NOT NULL default '0000-00-00 00:00:00',
  `logic_name` varchar(100) NOT NULL,
  `db` varchar(120) default NULL,
  `db_version` varchar(40) default NULL,
  `db_file` varchar(120) default NULL,
  `program` varchar(80) default NULL,
  `program_version` varchar(40) default NULL,
  `program_file` varchar(80) default NULL,
  `parameters` text,
  `module` varchar(80) default NULL,
  `module_version` varchar(40) default NULL,
  `gff_source` varchar(40) default NULL,
  `gff_feature` varchar(40) default NULL,
  PRIMARY KEY  (`analysis_id`),
  UNIQUE KEY `logic_name` (`logic_name`),
  KEY `logic_name_idx` (`logic_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `analysis_description`
--

DROP TABLE IF EXISTS `analysis_description`;
CREATE TABLE `analysis_description` (
  `analysis_id` smallint(5) unsigned NOT NULL,
  `description` text,
  `display_label` varchar(255) NOT NULL,
  `displayable` BOOLEAN NOT NULL default '1',
  `web_data` text,	
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `meta`
--

DROP TABLE IF EXISTS `meta`;
CREATE TABLE `meta` (
  `meta_id` int(10) NOT NULL auto_increment,
  `species_id` int(10) unsigned default '1',
  `meta_key` varchar(40) NOT NULL,
  `meta_value` varchar(950) NOT NULL,
  PRIMARY KEY  (`meta_id`),
  UNIQUE KEY `species_key_value_idx` (`species_id`,`meta_key`,`meta_value`),
  KEY `species_value_idx` (`species_id`,`meta_value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- Add empty schema_version entry for patches to update
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, "schema_version", "");
INSERT INTO meta (meta_key, meta_value) VALUES ('schema_type', 'funcgen');



--
-- Table structure for table `meta_coord`
--

DROP TABLE IF EXISTS `meta_coord`;
CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `max_length` int(11) default NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




--
-- Table structure for table `identity_xref`
--

DROP TABLE IF EXISTS identity_xref;
CREATE TABLE identity_xref (
  object_xref_id          INT(10) UNSIGNED NOT NULL,
  xref_identity 	  INT(5),	
  ensembl_identity        INT(5),
  xref_start              INT,
  xref_end                INT,
  ensembl_start           INT,
  ensembl_end             INT,
  cigar_line              TEXT, 
  score                   DOUBLE,
  evalue                  DOUBLE,
  PRIMARY KEY (object_xref_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `xref`
--

DROP TABLE IF EXISTS xref;
CREATE TABLE xref (
   xref_id 		      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
   external_db_id             SMALLINT UNSIGNED NOT NULL,
   dbprimary_acc              VARCHAR(40) NOT NULL,
   display_label              VARCHAR(128) NOT NULL,
   version                    VARCHAR(10) DEFAULT '0' NOT NULL,
   description                VARCHAR(255),
   info_type                  ENUM('PROJECTION', 'MISC', 'DEPENDENT', 'DIRECT', 'SEQUENCE_MATCH', 'INFERRED_PAIR', 'PROBE', 'UNMAPPED', 'CODING', 'TARGET') not NULL,
   info_text                  VARCHAR(255),
   PRIMARY KEY (xref_id),
   UNIQUE KEY id_index (dbprimary_acc, external_db_id, info_type, info_text),
   KEY display_index (display_label),
   KEY info_type_idx (info_type)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=100;


--
--  Table structure for table 'external_synonym'
--

DROP TABLE IF EXISTS external_synonym;
CREATE TABLE external_synonym (
  xref_id                     INT(10) UNSIGNED NOT NULL,
  synonym                     VARCHAR(40) NOT NULL, 
  PRIMARY KEY (xref_id, synonym),
  KEY name_index (synonym)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=20;


--
-- Table structure for table 'external_db' 
--

DROP TABLE IF EXISTS external_db;
CREATE TABLE external_db (
  external_db_id 	          SMALLINT(5) UNSIGNED NOT NULL auto_increment,
  db_name                     VARCHAR(100) NOT NULL,
  db_release                  VARCHAR(255),
  status                      ENUM('KNOWNXREF','KNOWN','XREF','PRED','ORTH', 'PSEUDO') NOT NULL,
  dbprimary_acc_linkable      BOOLEAN DEFAULT 1 NOT NULL,
  display_label_linkable      BOOLEAN DEFAULT 0 NOT NULL,
  priority                    INT NOT NULL,
  db_display_name             VARCHAR(255),
  type                        ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL') default NULL,
  secondary_db_name           VARCHAR(255) DEFAULT NULL,
  secondary_db_table          VARCHAR(255) DEFAULT NULL,
  description                 TEXT,
  PRIMARY KEY (external_db_id) 
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=80;


--
-- Table structure for table 'ontology_xref'
--

DROP TABLE if EXISTS ontology_xref;
CREATE TABLE ontology_xref (
  object_xref_id          INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  source_xref_id          INT(10) UNSIGNED DEFAULT NULL,
  linkage_type            ENUM('IC', 'IDA', 'IEA', 'IEP', 'IGI', 'IMP', 
		               'IPI', 'ISS', 'NAS', 'ND', 'TAS', 'NR', 'RCA') NOT NULL,
  KEY (object_xref_id),
  KEY (source_xref_id),
  UNIQUE (object_xref_id, source_xref_id, linkage_type)
)  ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- This is just an empty to table to avoid having to rework all the core sql and API to accomodate eFG specific xref schema

--
-- Table structure for table 'unmapped_reason'
--

DROP TABLE if EXISTS unmapped_reason;
CREATE TABLE `unmapped_reason` (
  `unmapped_reason_id` smallint(5) unsigned NOT NULL auto_increment,
  `summary_description` varchar(255) default NULL,
  `full_description` varchar(255) default NULL,
  PRIMARY KEY  (`unmapped_reason_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- CORE LIKE TABLES ---

--
-- Table structure for table `object_xref`
--

DROP TABLE IF EXISTS `object_xref`;
CREATE TABLE object_xref (
  object_xref_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  ensembl_id                  INT(10) UNSIGNED NOT NULL, 
  ensembl_object_type         ENUM('RegulatoryFeature', 'ExternalFeature', 'AnnotatedFeature', 'FeatureType', 'ProbeSet', 'Probe', 'ProbeFeature') not NULL,
  xref_id                     INT UNSIGNED NOT NULL,
  linkage_annotation          VARCHAR(255) DEFAULT NULL,
  analysis_id                 SMALLINT(5) UNSIGNED NOT NULL,
  UNIQUE (ensembl_object_type, ensembl_id, xref_id),
  KEY oxref_idx (object_xref_id, xref_id, ensembl_object_type, ensembl_id),
  KEY xref_idx (xref_id, ensembl_object_type),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=40;

-- Note we use case correct versions of object name to allow easy adaptor generation

--
-- Table structure for table `unmapped_object`
--

DROP TABLE IF EXISTS `unmapped_object`;
CREATE TABLE `unmapped_object` (
  `unmapped_object_id` int(10) unsigned NOT NULL auto_increment,
  `type` enum('xref', 'probe2transcript', 'array_mapping') NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `external_db_id` smallint(5) unsigned default NULL,
  `identifier` varchar(255) NOT NULL,
  `unmapped_reason_id` smallint(5) unsigned NOT NULL,
  `query_score` double default NULL,
  `target_score` double default NULL,
  `ensembl_id` int(10) unsigned default '0',
  `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType', 'Probe', 'ProbeSet', 'ProbeFeature') NOT NULL,
  `parent` varchar(255) default NULL,
  PRIMARY KEY  (`unmapped_object_id`),
  KEY `object_type_idx` (`ensembl_id`, `ensembl_object_type`),
  KEY `id_idx` (`identifier`),
  KEY `anal_idx` (`analysis_id`),
  KEY `anal_exdb_idx` (`analysis_id`,`external_db_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `coord_system`
--

DROP TABLE IF EXISTS `coord_system`;
CREATE TABLE `coord_system` (
  `coord_system_id` int(10) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL,
  `version` varchar(255) NOT NULL default '',
  `rank` int(11) NOT NULL,
  `attrib` set('default_version','sequence_level') default NULL,
  `schema_build` varchar(10) NOT NULL default '',
  `core_coord_system_id` int(10) NOT NULL,
  `species_id` int(10) NOT NULL default '1',
  `is_current` boolean default True,
  PRIMARY KEY  (`name`,`version`,`schema_build`,`species_id`),
  KEY `name_version_idx` (`name`,`version`),
  KEY `coord_species_idx` (`species_id`),
  KEY `coord_system_id_idx` (`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Could use boolean for is_current, but is MySQL implementation is non standard and an alias to tiny int anyway
-- This is only queried once to cache all the CSs on start up
-- Do we need attrib and rank for eFG?

--
-- Table structure for table `seq_region`
--

DROP TABLE IF EXISTS `seq_region`;
CREATE TABLE `seq_region` (
  `seq_region_id` int(10) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `core_seq_region_id` int(10) unsigned NOT NULL,
  `schema_build` varchar(10) NOT NULL default '',
  PRIMARY KEY  (`name`,`schema_build`,`coord_system_id`),
  KEY `coord_system_id` (`coord_system_id`),
  KEY `seq_region_id_idx` (`seq_region_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- it maybe possible to have 2 seq_regions on different levels with the same name
-- there fore have to use core cs id in primary key
-- order of extra primary key members doesn't really matter as we'll never query on them?
-- swapped order of ke to name coord_system_id, as we will probably want to query on name primarily
-- why does core table have cs id key?

-- Name is only required to enable us to add new seq_regions to the correct seq_region_id
-- It will never be used to retrieve a slice as we do that via the core DB
-- how are we going to use this when querying?
-- basically pull back seq_region_id based schema_build and core_seq_region_id
-- can we omit core_coord_system_id? As we have this info from the cs table.
-- other keys?





-- EFG tables




--
-- Table structure for table `associated_feature_type`
--



DROP TABLE IF EXISTS `associated_feature_type`;
CREATE TABLE `associated_feature_type` (
   `table_id` int(10) unsigned NOT NULL,
   `table_name` enum('annotated_feature', 'external_feature', 'regulatory_feature', 'feature_type') default NULL,
   `feature_type_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`table_id`, `table_name`, `feature_type_id`),
   KEY `feature_type_index` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




--
-- Table structure for table `experimental_group`
--

DROP TABLE IF EXISTS `experimental_group`;
CREATE TABLE `experimental_group` (
   `experimental_group_id` smallint(6) unsigned NOT NULL auto_increment,
   `name` varchar(40) NOT NULL,
   `location` varchar(120) default NULL,
   `contact` varchar(40) default NULL,
   PRIMARY KEY  (`experimental_group_id`),
   UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




--
-- Table structure for table `array`
--

DROP TABLE IF EXISTS `array`;
CREATE TABLE `array` (
   `array_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(40) default NULL,
   `format` varchar(20) default NULL,
   `vendor` varchar(40) default NULL,
   `description` varchar(255) default NULL,
   `type` varchar(20) default NULL,
   `class` varchar(20) default NULL,	
   PRIMARY KEY  (`array_id`),
   UNIQUE KEY  `vendor_name_idx` (`vendor`, `name`),
   UNIQUE KEY  `class_name_idx` (`class`, `name`)	
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- format = tiled, targetted, expression, custom/mixed?
-- class = AFFY_UTR AFFY_ST ILLUMINA_WG
-- Do we need to enum these?


--
-- Table structure for table `array_chip`
--

DROP TABLE IF EXISTS `array_chip`;
CREATE TABLE `array_chip` (
   `array_chip_id` int(10) unsigned NOT NULL auto_increment,
   `design_id` varchar(20) default NULL,
   `array_id` int(10) unsigned NOT NULL,
   `name` varchar(40) default NULL,
    PRIMARY KEY  (`array_chip_id`),
   UNIQUE KEY `array_design_idx` (`array_id`, `design_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- name = design_name
-- removed  `description` varchar(255) default NULL,


--
-- Table structure for table `probe_feature`
--

DROP TABLE IF EXISTS `probe_feature`;
CREATE TABLE `probe_feature` (
   `probe_feature_id` int(10) unsigned NOT NULL auto_increment,
   `seq_region_id` int(10) unsigned NOT NULL,
   `seq_region_start` int(10) NOT NULL,
   `seq_region_end` int(10) NOT NULL,
   `seq_region_strand` tinyint(4) NOT NULL, 
   `probe_id` int(10) unsigned NOT NULL,
   `analysis_id` smallint(5) unsigned NOT NULL,	
   `mismatches` tinyint(4) NOT NULL,
   `cigar_line` varchar(50) default NULL,
   PRIMARY KEY  (`probe_feature_id`),
   KEY `probe_idx` (`probe_id`),
   KEY `seq_region_probe_probe_feature_idx` (`seq_region_id`,`seq_region_start`, `seq_region_end`, `probe_id`, `probe_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `probe_set`
-- 

DROP TABLE IF EXISTS `probe_set`;
CREATE TABLE `probe_set` (
   `probe_set_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(100) NOT NULL,
   `size` smallint(6) unsigned NOT NULL,
   `family` varchar(20) default NULL,
   PRIMARY KEY  (`probe_set_id`),
	KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- family= ENCODE REGIONS, RANDOM etc, generic descriptor for probes, can have multiple families on one chip



--
-- Table structure for table `probe`
--

DROP TABLE IF EXISTS `probe`;
CREATE TABLE `probe` (
   `probe_id` int(10) unsigned NOT NULL auto_increment,
   `probe_set_id` int(10) unsigned default NULL,
   `name` varchar(40) NOT NULL,
   `length` smallint(6) unsigned NOT NULL,
   `array_chip_id` int(10) unsigned NOT NULL,
   `class` varchar(20) default NULL,
   `description` varchar(255) DEFAULT NULL,
    PRIMARY KEY  (`probe_id`, `name`, `array_chip_id`),
    KEY `probe_set_idx`  (`probe_set_id`),
    KEY `array_chip_idx` (`array_chip_id`),
	KEY `name_idx`       (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- remove? class = control, experimental etc... naming clash with array.class, different class types.
-- pair_index aka nimblegen match_index, id to connect paired probes (same as id/name for single probes).
-- 
-- Still have issues will true probe seq storage?  Some probes have mismatches and map multiple times therefore cannot use genome seq?
-- Seq storage not necessary for core, but maybe for FG analysis DB?
-- Have seq field/table?  Huge amount of probes and theoretically some could be >>1000bps.  If probes never overlapped would actually be less seq data as probe should be non-redundant within array.  Not the case, but may be easiest way to maintain true probe seq data.  Would negate the need for seq_region tables, this may hinder remapping. Then just have mismatches and cigar_line in probe_feature table instead of in both.  Would then also need build and chromosome tables.  How important is this?

-- As nimblegen probes are premapped to one location, can we capture this in probe_feature.analysis_id (e.g. vendor_mapping)?  Then alter this when we remap.  How can we maintain true match/mismatch pairs on remapping?  Only map PM probe, and then place MM probe at same location.  So mismatch in probe_feature would denote PM or MM probe only if probe/pair_index is valid and, vice versa.


-- Xref issue:  How are we going to consolidate vendor defined xrefs vs. ensembl core xrefs (e.g. affy) to ensure updating of xref table in core DB?



-- pair_index table   `pair_index` int(10) unsigned NOT NULL default '0',
-- joint index on probe_ids

-- Other fields:
--     length?
-- Removed:
--       `seq` mediumtext NOT NULL, - now captured via probe_feature.seq_region_id, see above for seq storage issues.
--	  `mismatch` tinyint(4) NOT NULL default '0', - Now using probe_feature, see above
--       array_id, now in probe_set


--
-- Table structure for table `probe_design`
--

DROP TABLE IF EXISTS `probe_design`;
CREATE TABLE `probe_design` (
   `probe_id` int(10) unsigned NOT NULL,
   `analysis_id` smallint(5) unsigned NOT NULL,
   `score` double default NULL,	
   `coord_system_id` int(10) unsigned NOT NULL,
    PRIMARY KEY  (`probe_id`, `analysis_id`, `coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- really wanted default NULL for coord_system_id but cannot have with int

--
-- Table structure for table `experiment`
--

DROP TABLE IF EXISTS `experiment`;
CREATE TABLE `experiment` (
   `experiment_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(100) default NULL,
   `experimental_group_id` smallint(6) unsigned default NULL,
   `date` date default '0000-00-00',
   `primary_design_type` varchar(30) default NULL, 
   `description`  varchar(255) default NULL,
   `mage_xml_id` int(10) unsigned default NULL,
   PRIMARY KEY  (`experiment_id`),
   UNIQUE KEY `name_idx` (`name`),
   KEY `design_idx` (`primary_design_type`),
   KEY `experimental_group_idx` (`experimental_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- alter table experiment add schema_build?  This would enable trace back to original import DB, just incse seq_region_ids change, or use status? IMPORTED_schema_build, MAPPED_schema_build or cs_id?
-- remove primary design_type
-- design_type  = CHIP2 etc... (is also design type in ontology i.e. binding_site_identification)
-- Secondary design type may be redundant, so have associated_design_types table, containing MGED ontology types?  Or have ontology_types table to control input and have linker table with just IDs?  Too normalised?


--
-- Table structure for table `experimental_design`
--

DROP TABLE IF EXISTS `experimental_design`;
CREATE TABLE `experimental_design` (
   `design_type_id` int(10) unsigned NOT NULL auto_increment,
   `table_name` varchar(40) default NULL,
   `table_id` int(10) unsigned default NULL,	
   PRIMARY KEY  (`design_type_id`, `table_name`, `table_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `mage_xml`
--

DROP TABLE IF EXISTS `mage_xml`;
CREATE TABLE `mage_xml` (
   `mage_xml_id` int(10) unsigned NOT NULL auto_increment,
   `xml` text,
   PRIMARY KEY  (`mage_xml_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;






-- Do we need a key just on table_id and name?

-- Handles design_types other than the primary experimental.primary_design_type

--
-- Table structure for table `design_type`
--

DROP TABLE IF EXISTS `design_type`;
CREATE TABLE `design_type` (
   `design_type_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(255) default NULL,
   PRIMARY KEY  (`design_type_id`),
   KEY `design_name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- add description here?
-- Handles design_types other than the primary experimental.primary_design_type



-- Table structure for table `feature_type`

DROP TABLE IF EXISTS `feature_type`;
CREATE TABLE `feature_type` (
   `feature_type_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(40) NOT NULL,
   `class` enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 'Regulatory Motif',  'Enhancer', 'Expression', 'Pseudo', 'Open Chromatin', 'Search Region', 'Association Locus') default NULL,
   `description`  varchar(255) default NULL,
   PRIMARY KEY  (`feature_type_id`),
   UNIQUE KEY `name_class_idx` (`name`, `class`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- Table to contain Brno nomenclature (modification ontology?) etc.
-- enum on class? HISTONE, PROMOTER
-- Have ontology name in field?
-- Have load ontology tool script?


-- Table structure for table `data_set`

DROP TABLE IF EXISTS `data_set`;
CREATE TABLE `data_set` (
   `data_set_id` int(10) unsigned NOT NULL auto_increment,
   `feature_set_id` int(10) unsigned default '0',
   `name` varchar(100) default NULL,
   PRIMARY KEY  (`data_set_id`, `feature_set_id`),
   UNIQUE KEY `name_idx` (name)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- PRIMARY KEY is just to impose uniqueness
-- any more keys?  Primary access would be through result_set_id to test if it appears as displayable in status, but this will just be a straight list of all
-- so maybe we can put that at the end of the index and the most used for query at the start?


-- Link table to provide many to one relationships:
--	target/feature < experiment (via experimental_chip_id)
--	prediction < experiment (via feature_set.feature_set_id)
-- Having this keyed experimental_chip also allows us to have a feature spanning two chips from the same exxperiment, aswell a delineating between potentially differentially formated chips from the same experiment.  Don't really need this level of info, just experiment_id will suffice.
-- Can we just use the index and drop the table?  Not that big so doesn't matter
-- Omitting experimental_feature_type_id from feature_set also allows us to drop this table for lite?


-- THis was originally also the record of the experiment feature_type, but this is now NR in this table.
-- We could move it to the result_set table, is more focused on the experimental_chips
-- But would potentially also be NR there too.  
-- Is still NR in result_set as we can have several analyses on the same chip/channels
-- Feature type is also known before results are imported on an exp/ec level, but would always import results anyway
-- So it's logically a little removed from where it should be, but more practical from a result set handling point of view.
-- could have in:
-- result_feature, but would not necessarily have feature_set_id
-- result_set
-- experimental_chip

-- we want to alter data_set to handle all type of data associations, rather than just feature < result
-- feature >< result would we ever get two feature sets from a group analysis?
-- feature < feature e.g. regulatory build/overlap sets
-- currently just have 1 to 1 relationship apart from reg build
-- let's just concentrate on later, many to many can be handle by duplicating data sets
-- ancilliary table with data set members, making data_set nr
-- can we have set type in data_set, rather than having redundant entries in data_set_member?
-- this will only work if we never have mixed types.




DROP TABLE IF EXISTS `supporting_set`;
CREATE TABLE `supporting_set` (
   `data_set_id` int(10) unsigned NOT NULL,
   `supporting_set_id` int(10) unsigned NOT NULL,
   `type` enum('result','feature','input') default NULL,
   PRIMARY KEY  (`data_set_id`, `supporting_set_id`),
   KEY `type_idx` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




-- Table structure for table `result`

DROP TABLE IF EXISTS `result`;
CREATE TABLE `result` (
   `result_id` int(10) unsigned NOT NULL auto_increment,
   `probe_id` int(10) unsigned default NULL,
   `score` double default NULL,
   `result_set_input_id` int(10) unsigned NOT NULL,
   `X` smallint(4) unsigned default NULL,
   `Y` smallint(4) unsigned default NULL,
   PRIMARY KEY  (`result_id`),
   KEY `probe_idx` (`probe_id`),
   KEY `result_set_input_idx` (`result_set_input_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

-- result_id needed as we may have replicate probe on same chip
-- X Y here allows repicate probes on same ship


-- Table structure for result_feature

DROP TABLE IF EXISTS `result_feature`;
CREATE TABLE `result_feature` (
  `result_feature_id` int(10) unsigned NOT NULL auto_increment,
  `result_set_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) NOT NULL,
  `seq_region_end` int(10) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `window_size` smallint(5) unsigned NOT NULL,
  `scores` longblob NOT NULL,
  KEY `result_feature_idx` (`result_feature_id`),
  KEY `set_window_seq_region_idx` (`result_set_id`, `window_size`,`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1
 PARTITION BY LIST (window_size)
 (PARTITION p0 VALUES IN (0),
 PARTITION p30 VALUES IN (30),
 PARTITION p65 VALUES IN (65),
 PARTITION p130 VALUES IN (130), 
 PARTITION p260 VALUES IN (260), 
 PARTITION p450 VALUES IN (450), 
 PARTITION p648 VALUES IN (648), 
 PARTITION p950 VALUES IN (950),
 PARTITION p1296 VALUES IN (1296)	
);

-- Partitions is set to number of windows
-- 0 partition optimised for natural res for arrays(0)
-- The rest are optimised for the default levels of zoom, plus a few intermediates 
-- No primary key as we need to parition on window_size, which would not be the whole of the primary key
-- set_window_region_idx can't be UNIQUE as there may be duplicates in the 0 wsize collections
-- i.e. two or more probe features with the same start, originating from replicate probes
-- or a probe seq which is a substr of another probe.

-- Table structure for `result_set`

DROP TABLE IF EXISTS `result_set`;
CREATE TABLE `result_set` (
   `result_set_id` int(10) unsigned NOT NULL auto_increment,
   `analysis_id` smallint(5) unsigned NOT NULL,
   `name` varchar(100) default NULL,
   `cell_type_id` int(10) unsigned default NULL,
   `feature_type_id` int(10) unsigned default NULL,
   PRIMARY KEY  (`result_set_id`),
   UNIQUE KEY `unique_idx` (`name`,`analysis_id`,`feature_type_id`,`cell_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `result_set_input`
--

DROP TABLE IF EXISTS `result_set_input`;
CREATE TABLE `result_set_input` (
   `result_set_input_id` int(10) unsigned NOT NULL auto_increment,
   `result_set_id` int(10) unsigned NOT NULL,
   `table_id` int(10) unsigned NOT NULL,
   `table_name` enum('experimental_chip','channel','input_set') DEFAULT NULL,
   PRIMARY KEY  (`result_set_input_id`, `result_set_id`),
   UNIQUE KEY `rset_table_idname_idx` (`result_set_id`, `table_id`, `table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




--
-- Table structure for table `annotated_feature`
--

DROP TABLE IF EXISTS `annotated_feature`;
CREATE TABLE `annotated_feature` (
  `annotated_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) default NULL,
  `score` double default NULL,
  `feature_set_id` int(10) unsigned NOT NULL,
  `summit` int(10) unsigned default NULL,
  PRIMARY KEY  (`annotated_feature_id`),
  UNIQUE KEY `seq_region_feature_set_idx` (`seq_region_id`,`seq_region_start`,`feature_set_id`),
  KEY `feature_set_idx` (`feature_set_id`)	  
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;



--
-- Table structure for table `motif_feature`
--

CREATE TABLE `motif_feature` (
  `motif_feature_id` int(10) unsigned NOT NULL auto_increment,
  `binding_matrix_id` INT(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) default NULL,
  `score` double default NULL,
  PRIMARY KEY  (`motif_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `binding_matrix_idx` (`binding_matrix_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `associated_motif_feature`
--

CREATE TABLE `associated_motif_feature` (
   `annotated_feature_id` int(10) unsigned NOT NULL,
   `motif_feature_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`annotated_feature_id`, `motif_feature_id`),
   KEY `motif_feature_idx` (`motif_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `binding_matrix`
--

CREATE  TABLE `binding_matrix` (
 `binding_matrix_id` INT(10) unsigned NOT NULL auto_increment,
 `name` VARCHAR(45) NOT NULL,
 `type` VARCHAR(45) NOT NULL,
 `feature_type_id` int(10) unsigned NOT NULL,
 `frequencies` VARCHAR(1000) NOT NULL,
 `description` VARCHAR(255) NULL,
 PRIMARY KEY (`binding_matrix_id`) ,
 KEY `name_type_idx` (`name`, `type`),
 KEY `feature_type_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- change frequencies to blob


-- Table structure for `feature_set`

DROP TABLE IF EXISTS `feature_set`;
CREATE TABLE `feature_set` (
   `feature_set_id` int(10) unsigned NOT NULL auto_increment,
   `feature_type_id` int(10) unsigned NOT NULL,
   `analysis_id` smallint(5) unsigned NOT NULL,
   `cell_type_id` int(10) unsigned default NULL,
   `name` varchar(100) default NULL,
   `type` enum('annotated', 'regulatory', 'external') default NULL,
   `description` varchar(80) default NULL,
   `display_label` varchar(80) default NULL,
   PRIMARY KEY  (`feature_set_id`),
   KEY `feature_type_idx` (`feature_type_id`),
   UNIQUE KEY `name_idx` (name)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `external_feature`
--

DROP TABLE IF EXISTS `external_feature`;
CREATE TABLE `external_feature` (
  `external_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,	
  `display_label` varchar(60) default NULL,
  `feature_type_id`	int(10) unsigned default NULL,
  `feature_set_id` int(10) unsigned NOT NULL,
  PRIMARY KEY  (`external_feature_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=80;




--
-- Table structure for table `regulatory_feature`
--


DROP TABLE IF EXISTS `regulatory_feature`;
CREATE TABLE `regulatory_feature` (
  `regulatory_feature_id` int(10) unsigned NOT NULL auto_increment,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(80) default NULL,
  `feature_type_id`	int(10) unsigned default NULL,
  `feature_set_id`	int(10) unsigned default NULL,
  `stable_id` mediumint(8) unsigned default NULL,
  `bound_seq_region_start` int(10) unsigned NOT NULL,	
  `bound_seq_region_end` int(10) unsigned NOT NULL,
  `binary_string` varchar(255) default NULL,
  `projected` boolean default FALSE,	
  PRIMARY KEY  (`regulatory_feature_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `stable_id_idx` (`stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;


-- Index on bound values?
-- stable_id is not unique, may have several instances across different cell_types/feature_sets


--
-- Table structure for table `regulatory_attribute`
--

DROP TABLE IF EXISTS `regulatory_attribute`;
CREATE TABLE `regulatory_attribute` (
  `regulatory_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated', 'external', 'motif') default NULL,
  PRIMARY KEY  (`regulatory_feature_id`, `attribute_feature_table`, `attribute_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;


-- Cell specific build:
-- Would need to add cell_type_id here? Or maybe we could simply have another feature_set representing
-- each cell_type specific build.
-- The feature_set_id in regulatory_feature table would be the main 'RegulatoryFeatures' feature_set
-- and the feature_set_ids in regulatory_attribute would be e.g. 'GM06690 - RegulatoryFeatures'
-- This would not behave like a normal feature_set, as we would not be able to retrieve it like the others
-- Maybe we could give it a different feature_set type, and code the FeatureSetAdaptor/RegulatoryFeatureAdaptor to handle this
-- Or do we need a slightly different class?



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

DROP TABLE IF EXISTS `input_set`;
CREATE TABLE `input_set` (
   `input_set_id` int(10) unsigned NOT NULL auto_increment,
   `experiment_id` int(10) unsigned default NULL,
   `feature_type_id` int(10) unsigned default NULL,
   `cell_type_id` int(10) unsigned default NULL,
   `format` varchar(20) default NULL,
   `vendor` varchar(40) default NULL,
   `name` varchar(100) not NULL,
   `type` enum('annotated', 'result') default NULL,
   PRIMARY KEY  (`input_set_id`),
   UNIQUE KEY `name_idx` (`name`),
   KEY `experiment_idx` (`experiment_id`),
   KEY `feature_type_idx` (`feature_type_id`),
   KEY `cell_type_idx` (`cell_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=30;

-- keys on vendor/format/type?
-- type = annotated, result
-- What about defining coverage? To determine if we don't want to load 0 window_size
-- for high volume seq align sets? Just do this during import
-- But we need to know whether a set has 0 or not for retreiving the ResultFeatures
-- Do we have window_size table?
-- Could also potentially have float scores in Input, so need to capture that too for unpacking?


--
-- Table structure for table `input_subset`
--

-- represents a file from a subset
-- mainly used for tracking import and recovery 

DROP TABLE IF EXISTS `input_subset`;
CREATE TABLE `input_subset` (
   `input_subset_id` int(10) unsigned NOT NULL auto_increment,
   `input_set_id` int(10) unsigned NOT NULL,
   `name` varchar(100) NOT NULL, -- filename?	
   PRIMARY KEY  (`input_subset_id`), 
   UNIQUE KEY `set_name_dx` (`input_set_id`, `name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=30;


-- remove experimental_set_id from key to force uniqueness of file?
-- No file names may be same across expeirments, this allows addition of same sub set to different sets
-- whilst making it unique whtin a set
-- no we can't have multi experiment sets at present, still migh have same filename tho

--
-- Table structure for table `experimental_chip`
--

DROP TABLE IF EXISTS `experimental_chip`;
CREATE TABLE `experimental_chip` (
   `experimental_chip_id` int(10) unsigned NOT NULL auto_increment,
   `unique_id` varchar(20) NOT NULL,
   `experiment_id` int(10) unsigned default NULL,
   `array_chip_id` int(10) unsigned default NULL,
   `feature_type_id` int(10) unsigned default NULL,
   `cell_type_id` int(10) unsigned default NULL,
   `biological_replicate` varchar(100) default NULL,
   `technical_replicate` varchar(100) default NULL,
   PRIMARY KEY  (`experimental_chip_id`),
   KEY `experiment_idx` (`experiment_id`),
   KEY `feature_type_idx` (`feature_type_id`),
   KEY `unique_id_idx` (`unique_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- We can't implement uniqueness of ec(unique_id for a given vendor) via key here
-- There may be clashes between vendors, and a chip could potentially be re-used in another exp
-- This is handled by the fetch methods:
--    fetch_by_unique_and_experiment_id
--    fetch_by_unique_id_vendor
-- add key on cell_type? or remove key on feature_type.  Will we ever want to query based on these? I doubt it.
-- Should handle re-usage of physical chip


--
-- Table structure for table `channel`
--

DROP TABLE IF EXISTS `channel`;
CREATE TABLE `channel` (
   `channel_id` int(10) unsigned NOT NULL auto_increment,
   `experimental_chip_id` int(10) unsigned default NULL,	
   `sample_id` varchar(20) default NULL,
   `dye`  varchar(20) default NULL,
   `type` varchar(20) default NULL,
   PRIMARY KEY  (`channel_id`),
   KEY `experimental_chip_idx` (`experimental_chip_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- type should be restricted to EXPERIMENTAL & CONTROL?
-- all other variables should be in experimental_variable table

--
-- Table structure for table `cell_type`
--

DROP TABLE IF EXISTS `cell_type`;
CREATE TABLE `cell_type` (
   `cell_type_id` int(10) unsigned NOT NULL auto_increment,
   `name`  varchar(120) not NULL,
   `display_label` varchar(20) default NULL,
   `description` varchar(80) default NULL,
   `gender` enum('male', 'female') default NULL,
   PRIMARY KEY  (`cell_type_id`),
   UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- type? enum Tissue? Line?
-- xref to coriell?

--
-- Table structure for table `experimental_variable`
--

DROP TABLE IF EXISTS `experimental_variable`;
CREATE TABLE `experimental_variable` (
   `table_id` int(10) unsigned default NULL,
   `table_name` varchar(20) default NULL,	
   `name` varchar(40) default NULL,
   `unit` varchar(40) default NULL,
   `value` varchar(40) default NULL,
   PRIMARY KEY  (`table_id`, `table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- context_id = table_name:table_name_id
-- Needs to accommodate QC flag and description, and reusing chips.
-- Have list of mandatory defs for each group, then also test MGED ontology and warn if not present.
-- status for chip?
-- Some of these are chip rather than channels specific
-- `tissue` varchar(40) NOT NULL default '',
-- `cell_type` varchar(40) NOT NULL default '',
-- `dev_stage` varchar(40) NOT NULL default '',
-- `antibody` varchar(40) NOT NULL default '',
-- `time_point` time NOT NULL default '000:00:00',
-- antibody?
-- species
-- description?
-- can this handle experimental management such that slides can be reused and classed as such
-- Or classed as failures etc.
-- Do we also need an extra table to hold experiment level meta data along side design type, and rename this channel variable

--
-- Table structure for table `status`
--

DROP TABLE IF EXISTS `status`;
CREATE TABLE `status` (
   `table_id` int(10) unsigned default NULL,
   `table_name` varchar(32) default NULL,	
   `status_name_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`table_id`, `table_name`, `status_name_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- 32 is max length for mysql table name

--
-- Table structure for table `status`
--

DROP TABLE IF EXISTS `status_name`;
CREATE TABLE `status_name` (
   `status_name_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(20) default NULL,	
   PRIMARY KEY  (`status_name_id`),
   UNIQUE KEY `status_name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



-- Remove these to separate file and handle with import_type.pl?
INSERT into status_name(name) values ('DISPLAYABLE');
INSERT into status_name(name) values ('IMPORTED');
INSERT into status_name(name) values ('DAS_DISPLAYABLE');
INSERT into status_name(name) values ('RESOLVED');
INSERT into status_name(name) values ('VSN_GLOG');
INSERT into status_name(name) values ('Parzen');
INSERT into status_name(name) values ('T.Biweight');
INSERT into status_name(name) values ('LOESS');
INSERT into status_name(name) values ('MART_DISPLAYABLE');
INSERT into status_name(name) values ('RESULT_FEATURE_SET');
-- These are now species specific and should be added like the feature/cell types 
INSERT into status_name(name) values ('IMPORTED_NCBI36');
INSERT into status_name(name) values ('IMPORTED_GRCh37');


-- need to add more states, probably need to validate/insert required states in Importer
-- would need to get CoordSys objects and set IMPORTED_CS_"cs_id" for relevant data_version




-- Further thoughts:


-- Denormalise --
-- To mart style to optimize queries? Mart style interface to export(to R)?
-- probe_set, probe, probe_feature and results could all e split on array vendor, format or experiment.design_type
--     probe > "array.vendor"_probe e.g. affy_probe, or "array.name"_probe e.g u133_probe, or "array.class"_probe e.g. CHIP2_probe


