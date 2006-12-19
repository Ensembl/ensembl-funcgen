--mysql -hia64g -uensadmin -pensembl < efg_test.sql

--- Summary of current proposed/absent tables:
---     seq_region tables (incorporating chromosome name and assembly version)
---     xref tables - for linking probes to user/vendor defined features.
---     Denormalised probe_set/probe/probe_feature/results tables bt vendor, format or design_type? 
---     egroup_member table? 


--- Mapping between assemblies?
--- All default values need checking
--- Use MGED terms for design_types & experimental_variables where possible.
--- Use Brno nomelcature for target features



--
-- Table structure for table `egroup`
--

DROP TABLE IF EXISTS `egroup`;
CREATE TABLE `egroup` (
   `egroup_id` smallint(6) unsigned NOT NULL auto_increment,
   `name` varchar(40) NOT NULL default '',
   `location` varchar(120) default NULL,
   `contact` varchar(40) default NULL,
   PRIMARY KEY  (`egroup_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


#insert into egroup values("", "efg", "Hinxton", "njohnson@ebi.ac.uk");

--- group is reserved by MySQL.
---Others? head? description? egroup_member(table?) probably overkill.


--- Need to separate fields between array and array_chip sensibly
--- It seems Nimblgen can do pretty much anything they like in a chip set, but can we make any safe assumption that will be array/chip set specific

--
-- Table structure for table `array`
--

DROP TABLE IF EXISTS `array`;
CREATE TABLE `array` (
   `array_id` int(11) unsigned NOT NULL auto_increment,
   `name` varchar(40) default NULL,
   `format` varchar(20) default NULL,
   `vendor` varchar(40) default NULL,
   `description` varchar(255) default NULL,
   PRIMARY KEY  (`array_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- removed `size` tinyint(4) unsigned NOT NULL default '0', now dynamically generated from key array_chips
--- Size -> Chips?  Need some level of validation for chip sets to avoid incomplete chip sets
--- name = design_name, or is this the chip name?
--- format = tiled, gene, exon, targetted, custom/mixed? Need control/restrict these etc..
--- species, could we have multi-species arrays?
--- removed:
--   `design_id` int(11) unsigned NOT NULL default '0',

--
-- Table structure for table `array_chip`
--

DROP TABLE IF EXISTS `array_chip`;
CREATE TABLE `array_chip` (
   `array_chip_id` int(11) unsigned NOT NULL auto_increment,
   `design_id` varchar(20) default NULL,
   `array_id` int(11) unsigned NOT NULL,
   `name` varchar(40) default NULL,
    PRIMARY KEY  (`array_chip_id`),
   KEY `design_idx` (`design_id`),
   KEY `array_idx` (`array_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



--- name = design_name, or is this the chip name?
--- removed  `description` varchar(255) default NULL,




--
-- Table structure for table `oligo_feature`
--

DROP TABLE IF EXISTS `oligo_feature`;
CREATE TABLE `oligo_feature` (
   `oligo_feature_id` int(11) unsigned NOT NULL auto_increment,
   `seq_region_id` int(11) unsigned NOT NULL default '0',
   `seq_region_start` int(11) NOT NULL default '0',
   `seq_region_end` int(11) NOT NULL default '0',
   `seq_region_strand` tinyint(4) NOT NULL default '0', 
   `coord_system_id` int(10) unsigned NOT NULL default '0',
   `oligo_probe_id` int(11) unsigned NOT NULL default '0',
   `analysis_id` int(11) unsigned NOT NULL default '0',	
   `mismatches` tinyint(4) NOT NULL default '0',
   `cigar_line` text,
   PRIMARY KEY  (`oligo_feature_id`),
   KEY `oligo_probe_idx` (`oligo_probe_id`),
   KEY `seq_region_idx` (`seq_region_id`, `seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- remove mismatches?
--- Do jointindex on seq region local start? ADD COORD_SYSTEM to joint index?
--- Currently use chr data in here rather than true seq_regions, Need to map all probe global values to seq_regions.
--- build_id currently set to freeze date, can go when we incorporate build mapping into the import API
--- Would not capture true probe seq if there were any mismatches in mapping.
--- Mismatches in mapping, should default be NULL?  Mismatch for nimblegen here would denote MM of PM/MM probe pair
--- Cigarline(from alignment) shows mismatches, would normally be "probe.length"M, would need original probe seq for alignment, not currently stored with seq_region style tables.  See comments below probe table re: seq storage.
--- Joint index on seq_region_id and local start? 



--
-- Table structure for table `oligo_probe_set`
-- 

DROP TABLE IF EXISTS `oligo_probe_set`;
CREATE TABLE `oligo_probe_set` (
   `oligo_probe_set_id` int(11) unsigned NOT NULL auto_increment,
   `name` varchar(20) NOT NULL default '',
   `size` smallint(6) unsigned NOT NULL default '0',
   `family` varchar(20) default NULL,
   PRIMARY KEY  (`oligo_probe_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--- Now ancillary/optional table
--- joint key on probe_set_id and array_id?
--- aka feature for nimblegen(optional, therefore some probes = their featureset/probeset)
--- Can we omit probe_sets for sets of size == 1?
--- family= ENCODE REGIONS, RANDOM etc, generic descriptor for probes, can have multiple families on one chip
--- xref_id = encode region name, removed!!!!  Generate real xref entries.


---
-- Table structure for table `oligo_probe`
--

DROP TABLE IF EXISTS `oligo_probe`;
CREATE TABLE `oligo_probe` (
   `oligo_probe_id` int(11) unsigned NOT NULL auto_increment,
   `oligo_probe_set_id` int(11) unsigned default NULL,
   `name` varchar(40) NOT NULL default '',
   `length` smallint(6) unsigned NOT NULL default '0',
   `array_chip_id` int(11) unsigned NOT NULL default '0',
   `class` varchar(20) default NULL,
    PRIMARY KEY  (`oligo_probe_id`, `name`),
    KEY `probe_set_idx` (`oligo_probe_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--- remove array_id and link through probe_set, can't if we're not populating probe_set for sets of 1.
--- name (for Affy) at least is array_chip to probe relationship, probeset is optional for some formats 
--- remove? class = control, experimental etc... naming clash with array.class, different class types.
--- pair_index aka nimblegen match_index, id to connect paired probes (same as id/name for single probes).
--- 
--- Still have issues will true probe seq storage?  Some probes have mismatches and map multiple times therefore cannot use genome seq?
--- Seq storage not necessary for core, but maybe for FG analysis DB?
--- Have seq field/table?  Huge amount of probes and theoretically some could be >>1000bps.  If probes never overlapped would actually be less seq data as probe should be non-redundant within array.  Not the case, but may be easiest way to maintain true probe seq data.  Would negate the need for seq_region tables, this may hinder remapping. Then just have mismatches and cigar_line in probe_feature table instead of in both.  Would then also need build and chromosome tables.  How important is this?

--- As nimblegen probes are premapped to one location, can we capture this in probe_feature.analysis_id (e.g. vendor_mapping)?  Then alter this when we remap.  How can we maintain true match/mismatch pairs on remapping?  Only map PM probe, and then place MM probe at same location.  So mismatch in probe_feature would denote PM or MM probe only if probe/pair_index is valid and, vice versa.


--- Xref issue:  How are we going to consolidate vendor defined xrefs vs. ensembl core xrefs (e.g. affy) to ensure updating of xref table in core DB?



--- pair_index table   `pair_index` int(11) unsigned NOT NULL default '0',
--- joint index on oligo_probe_ids




--- Other fields:
---     length?
--- Removed:
---       `seq` mediumtext NOT NULL, - now captured via probe_feature.seq_region_id, see above for seq storage issues.
---	  `mismatch` tinyint(4) NOT NULL default '0', - Now using probe_feature, see above
---       array_id, now in probe_set
--
-- Table structure for table `experiment`
--

DROP TABLE IF EXISTS `experiment`;
CREATE TABLE `experiment` (
   `experiment_id` int(11) unsigned NOT NULL auto_increment,
   `name` varchar(30) default NULL,
   `egroup_id` smallint(6) unsigned default NULL,
   `date` date NOT NULL default '0000-00-00',
   `primary_design_type` varchar(30) default NULL, 
   `description`  varchar(255) default NULL,
   PRIMARY KEY  (`experiment_id`),
   KEY `egroup_idx` (`egroup_id`),
   KEY `design_idx` (`primary_design_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- joint key on egroup and array?
-- remove primary design_type
--- design_type  = CHIP2 etc... (is also design type in ontology i.e. binding_site_identification)
--- Secondary design type may be redundant, so have associated_design_types table, containing MGED ontology types?  Or have ontology_types table to control input and have linker table with just IDs?  Too normalised?
--- epi_feature would also be redundant for non-CHIP2 experiments...have associated_target table?  Or just have target field, with "Experssion" for standard Xn chips?


--- Other fields:
---     egroup_member? (overkill?)

--- removed:
---    `array_id` int(11) unsigned default NULL, Cannot guarantee one array/chip set, especially as we may have to capture Nimblegen chips as individual chip sets.


--
-- Table structure for table `experimental_design`
--

DROP TABLE IF EXISTS `experimental_design`;
CREATE TABLE `design_type` (
   `design_type_id` int(11) unsigned NOT NULL auto_increment,
   `table_name` varchar(40) default NULL,
   `table_id` int(11) unsigned default NULL,	
   PRIMARY KEY  (`design_type_id`, `table_name`, `table_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



-- Do we need a key just on table_id and name?

-- Handles design_types other than the primary experimental.primary_design_type

--
-- Table structure for table `design_type`
--

DROP TABLE IF EXISTS `design_type`;
CREATE TABLE `design_type` (
   `design_type_id` int(11) unsigned NOT NULL auto_increment,
   `name` varchar(255) default NULL,
   PRIMARY KEY  (`design_type_id`),
   KEY `design_name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- add description here?
-- Handles design_types other than the primary experimental.primary_design_type



-- Table structure for table `feature_type`

DROP TABLE IF EXISTS `feature_type`;
CREATE TABLE `feature_type` (
   `feature_type_id` int(11) unsigned NOT NULL auto_increment,
   `name` varchar(40) default NULL,
   `class` varchar(40) default NULL,
   `description`  varchar(255) default NULL,
   PRIMARY KEY  (`feature_type_id`),
   KEY `feature_type_name_class_idx` (`name`, `class`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--- Table to contain Brno nomenclature (modification ontology?) etc.
--- enum on class? HISTONE, PROMOTER
--- Have ontology name in field?
--- Have load ontology tool script?


-- Table structure for table `result_feature`

DROP TABLE IF EXISTS `result_feature`;
CREATE TABLE `result_feature` (
   `result_set_id` int(11) unsigned NOT NULL,
   `feature_set_id` int(11) unsigned NOT NULL,
   PRIMARY KEY  (`result_set_id`, `feature_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Link table to provide many to one relationships:
--	target/feature < experiment (via experimental_chip_id)
--	prediction < experiment (via feature_set.feature_set_id)
-- Having this keyed experimental_chip also allows us to have a feature spanning two chips from the same exxperiment, aswell a delineating between potentially differentially formated chips from the same experiment.  Don't really need this level of info, just experiment_id will suffice.
-- Can we just use the index and drop the table?  Not that big so doesn't matter
-- Omitting experimental_feature_type_id from feature_set also allows us to drop this table for lite?


--THis was originally also the record of the experiment feature_type, but this is now NR in this table.
-- We could move it to the result_set table, is more focused on the experimental_chips
-- But would potentially also be NR there too.  
-- Is still NR in result_set as we can have several analyses on the same chip/channels
-- Feature type is also known before results are imported on an exp/ec level, but would always import results anyway
-- So it's logically a little removed from where it should be, but more practical from a result set handling point of view.
-- could have in:
---- result_feature, but would not necessarily have feature_set_id
---- result_set
---- experimental_chip


-- Table structure for table `experiment_prediction`

--DROP TABLE IF EXISTS `experiment_prediction`;
--CREATE TABLE `experiment_prediction` (
--  `experiment_id` int(11) unsigned default NULL,
--   `predicted_feature_id` int(11) unsigned default NULL,
--   PRIMARY KEY  (`experiment_id`, `predicted_feature_id`)
--) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Superceded by feature_set
-- Provides link between predicted_features and contributing experiments
-- Which may be a subset of those in experiment_target

-- Table structure for table `result`

DROP TABLE IF EXISTS `result`;
CREATE TABLE `result` (
   `result_id` int(11) unsigned NOT NULL auto_increment,
   `oligo_probe_id` int(11) unsigned default NULL,
   `score` double default NULL,
   `analysis_id` int(11) unsigned default NULL,
   `table_id` int(11) unsigned default NULL,
   `table_name` varchar(20) default NULL,
   PRIMARY KEY  (`result_id`),
   KEY `oligo_probe_idx` (`oligo_probe_id`),
   KEY `table_name_id_idx` (`table_name`, `table_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000  AVG_ROW_LENGTH=40;


--- joint index on all but result_id and maybe analysis_id?
--- REMOVEd experimental_id?   `experimental_chip_id` int(11) unsigned default NULL,
--- joint primary key with probe_feature_id?
--- metric default would be id for "RAW"..no, need to test has been specifically set, so NULL
--- Allows storage of none raw values
---Also needs to accommodate different normalisations 


--- Table structure for `result_set`

DROP TABLE IF EXISTS `result_set`;
CREATE TABLE `result_set` (
   `result_set_id` int(11) unsigned NOT NULL auto_increment,
   `analysis_id` int(11) unsigned default NULL,
   `table_id` int(11) unsigned default NULL,
   `table_name` varchar(20) default NULL,
   `chip_set_id` smallint(6) unsigned default NULL,
   PRIMARY KEY  (`result_set_id`),
   KEY `table_name_id_idx` (`table_name`, `table_id`),
   KEY `chip_set_idx` (`chip_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--- Need to implement this




--
-- Table structure for table `predicted_feature`
--

DROP TABLE IF EXISTS `predicted_feature`;
CREATE TABLE `predicted_feature` (
  `predicted_feature_id` int(11) unsigned NOT NULL auto_increment,
  `seq_region_id` int(11) unsigned NOT NULL default '0',
  `seq_region_start` int(11) unsigned NOT NULL default '0',
  `seq_region_end` int(11) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(1) NOT NULL default '0',
  `coord_system_id` int(11) unsigned NOT NULL default '0',
  `feature_type_id` int(11) unsigned NOT NULL default '0',
  `feature_set_id` int(11) unsigned NOT NULL default '0',	
  `display_label` varchar(60) NOT NULL default '',
  `analysis_id` int(11) unsigned NOT NULL default '0',
  `score` double default NULL,
  PRIMARY KEY  (`predicted_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `type_idx` (`feature_type_id`),	  
  KEY `hit_idx` (`display_label`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

--- Need to be able to maintain link between prediction and source data i.e. experiment/s
--- Need to be able to have this table running completely from the persistent table
--- i.e. We want to be able to remove the experiment table data and still have access to all relevant data
--- e.g. primary design type, target_name/description 
--- Either make experiment->target 1:1 and have mixed persistent data in target(rm experiment_id and add target_id to experment)
--- Or add target_id to predicted_feature, may be redundant for non-tiling arrays?  How are we going to capture primary_design_type?


--- Table structure for `feature_set`

DROP TABLE IF EXISTS `feature_set`;
CREATE TABLE `feature_set` (
   `feature_set_id` int(11) unsigned NOT NULL auto_increment,
   `feature_type_id` int(11) unsigned NOT NULL,
   `analysis_id`  int(11) unsigned default NULL,
   `cell_type_id` int(11) unsigned default NULL,
   PRIMARY KEY  (`feature_set_id`),
   KEY `feature_type_idx` (`feature_type_id`)	
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- remove displayable
--- Need to implement this
-- remove displayable index?
-- change default NOT/NULLs?
-- feature_type_id is NR between this and pf, but need it here for ft type ResultSet queries 



--
-- Table structure for table `experimental_chip`
--

DROP TABLE IF EXISTS `experimental_chip`;
CREATE TABLE `experimental_chip` (
   `experimental_chip_id` int(11) unsigned NOT NULL auto_increment,
   `unique_id` varchar(20) NOT NULL default '0',
   `experiment_id` int(11) unsigned default NULL,
   `array_chip_id` int(11) unsigned default NULL,
   `feature_type_id` int(11) unsigned default NULL,
   `cell_type_id` int(11) unsigned default NULL,
   `description` varchar(255) default NULL,
   PRIMARY KEY  (`experimental_chip_id`),
   KEY `experiment_idx` (`experiment_id`),
   KEY `feature_type_idx` (`feature_type_id`),
   KEY `chip_idx` (`unique_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- add cell line key?
-- composite unique key onchi/exp/arra_chip ids?
--Should handle re-usage of physical chip
--Rename slide? or have array_chip, and experimental_chip

--
-- Table structure for table `channel`
--

DROP TABLE IF EXISTS `channel`;
CREATE TABLE `channel` (
   `channel_id` int(11) unsigned NOT NULL auto_increment,
   `experimental_chip_id` int(11) unsigned default NULL,	
   `sample_id` varchar(20) default NULL,
   `dye`  varchar(20) default NULL,
   `type` varchar(20) default NULL,
   `description` varchar(255) default NULL,
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
   `cell_type_id` int(11) unsigned NOT NULL auto_increment,
   `name`  varchar(120) default NULL,
   `display_label` varchar(20) default NULL,
   PRIMARY KEY  (`cell_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



--
-- Table structure for table `experimental_variable`
--

DROP TABLE IF EXISTS `experimental_variable`;
CREATE TABLE `experimental_variable` (
   `table_id` int(11) unsigned default NULL,
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
--- can this handle experimental management such that slides can be reused and classed as such
--- Or classed as failures etc.
--- Do we also need an extra table to hold experiment level meta data along side design type, and rename this channel variable

--
-- Table structure for table `status`
--

DROP TABLE IF EXISTS `status`;
CREATE TABLE `status` (
   `table_id` int(11) unsigned default NULL,
   `table_name` varchar(20) default NULL,	
   `state` varchar(40) default NULL,
   PRIMARY KEY  (`table_id`, `table_name`, `state`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



--
-- Table structure for table `analysis`
--

DROP TABLE IF EXISTS `analysis`;
CREATE TABLE `analysis` (
  `analysis_id` int(10) unsigned NOT NULL auto_increment,
  `created` datetime NOT NULL default '0000-00-00 00:00:00',
  `logic_name` varchar(40) NOT NULL default '',
  `db` varchar(120) default NULL,
  `db_version` varchar(40) default NULL,
  `db_file` varchar(120) default NULL,
  `program` varchar(80) default NULL,
  `program_version` varchar(40) default NULL,
  `program_file` varchar(80) default NULL,
  `parameters` varchar(255) default NULL,
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
  `analysis_id` int(10) unsigned NOT NULL default '0',
  `description` text,
  `display_label` varchar(255) default NULL,
  `displayable` tinyint(1) NOT NULL default '1',
  `web_data` text,	
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `meta`
--

DROP TABLE IF EXISTS `meta`;
CREATE TABLE `meta` (
  `meta_id` int(11) NOT NULL auto_increment,
  `meta_key` varchar(40) NOT NULL default '',
  `meta_value` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`meta_id`),
  KEY `meta_key_index` (`meta_key`),
  KEY `meta_value_index` (`meta_value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--
-- Table structure for table `meta_coord`
--

DROP TABLE IF EXISTS `meta_coord`;
CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL default '',
  `coord_system_id` int(11) NOT NULL default '0',
  `max_length` int(11) default NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--- change to primary key?
---should only ever be predicted_feature, but with all the coord_sys_ids
---This is slightly redundant, but required for core API modules to work

--- Set up default meta coord entries...this should be done in import
-- Change max lenght of oligo?

insert into meta_coord values("predicted_feature", 1, 147);
insert into meta_coord values("oligo_feature", 1, 50);

--
-- Table structure for table `coord_system`
--

DROP TABLE IF EXISTS `coord_system`;
CREATE TABLE `coord_system` (
  `coord_system_id` int(11) NOT NULL auto_increment,
  `name` varchar(40) NOT NULL default '',
  `version` varchar(40) default NULL,
  `rank` int(11) NOT NULL default '0',
  `attrib` set('default_version','sequence_level') default NULL,
  `schema_build` varchar(6) default NULL,
  PRIMARY KEY  (`coord_system_id`),
  UNIQUE KEY `rank` (`rank`, `schema_build`),
  UNIQUE KEY `name` (`name`,`version`, `schema_build`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

---Should only ever be chromosome?

--- Further thoughts:

--- More tables required for probe remapping? meta, rule tables etc?


--- Denormalise ---
--- To mart style to optimize queries? Mart style interface to export(to R)?
--- probe_set, probe, probe_feature and results could all e split on array vendor, format or experiment.design_type
---     probe > "array.vendor"_probe e.g. affy_probe, or "array.name"_probe e.g u133_probe, or "array.class"_probe e.g. CHIP2_probe


--- Mostly empty/unused fields which could be extracted to separate tables to reduce size?
---     cigar_line
