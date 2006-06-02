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

DROP DATABASE IF EXISTS `efg_test`;
CREATE DATABASE `efg_test`;
USE `efg_test`;

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
   `size` int(11) unsigned NOT NULL default '0',
   `species` varchar(20) default NULL,
   `vendor` varchar(40) default NULL,
   `description` varchar(255) default NULL,
   PRIMARY KEY  (`array_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


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
   `design_id` int(11) unsigned NOT NULL default '0',
   `array_id` int(11) unsigned NOT NULL default '0',
   `name` varchar(40) default NULL,
   `description` varchar(255) default NULL,
   PRIMARY KEY  (`array_chip_id`),
   KEY `design_idx` (`design_id`),
   KEY `array_idx` (`array_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



--- name = design_name, or is this the chip name?




--
-- Table structure for table `probe_feature`
--

DROP TABLE IF EXISTS `probe_feature`;
CREATE TABLE `probe_feature` (
   `probe_feature_id` int(11) unsigned NOT NULL auto_increment,
   `seq_region_id` int(11) unsigned NOT NULL default '0',
   `local_start` int(11) NOT NULL default '0',
   `local_end` int(11) NOT NULL default '0',
   `seq_region_strand` tinyint(4) NOT NULL default '0', 
   `mismatches` tinyint(4) NOT NULL default '0',
   `probe_id` int(11) unsigned NOT NULL default '0',
   `build_id` smallint(6) unsigned NOT NULL default '0',
   `analysis_id` int(11) unsigned NOT NULL default '0',
   `cigar_line` text,
   PRIMARY KEY  (`probe_feature_id`),
   KEY `start_idx` (`local_start`),
   KEY `probe_idx` (`probe_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--- Currently use chr data in here rather than true seq_regions, Need to map all probe global values to seq_regions.
--- build_id currently set to freeze date, can go when we incorporate build mapping into the import API
--- Would not capture true probe seq if there were any mismatches in mapping.
--- Mismatches in mapping, should default be NULL?  Mismatch for nimblegen here would denote MM of PM/MM probe pair
--- Cigarline(from alignment) shows mismatches, would normally be "probe.length"M, would need original probe seq for alignment, not currently stored with seq_region style tables.  See comments below probe table re: seq storage.
--- Joint index on seq_region_id and local start? 



--
-- Table structure for table `probe_set`
-- 

DROP TABLE IF EXISTS `probe_set`;
CREATE TABLE `probe_set` (
   `probe_set_id` int(11) unsigned NOT NULL auto_increment,
   `name` varchar(20) NOT NULL default '',
   `size` smallint(6) unsigned NOT NULL default '0',
   `array_chip_id` int(11) unsigned NOT NULL default '0',
   `family` varchar(20) default NULL,
   `xref_id` int(10) unsigned NOT NULL default '0',
   PRIMARY KEY  (`probe_set_id`),
   KEY `array_chip_idx` (`array_chip_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--- joint key on probe_set_id and array_id?
--- aka feature for nimblegen(optional, therefore some probes = their featureset/probeset)
--- Can we omit probe_sets for sets of size == 1?
--- family= ENCODE REGIONS, RANDOM etc, generic descriptor for probes, can have multiple families on one chip
--- xref_id = encode region name
-- Table structure for table `probe`
--

DROP TABLE IF EXISTS `probe`;
CREATE TABLE `probe` (
   `probe_id` int(11) unsigned NOT NULL auto_increment,
   `probe_set_id` int(11) unsigned NOT NULL default '0',
   `name` varchar(20) default NULL,
   `pair_index` int(11) unsigned NOT NULL default '0',
   `class` varchar(20) default NULL,
    PRIMARY KEY  (`probe_id`),
    KEY `probe_set_idx` (`probe_set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--- remove array_id and link through probe_set, can't if we're not populating probe_set for sets of 1.
--- remove? class = control, experimental etc... naming clash with array.class, different class types.
--- pair_index aka nimblegen match_index, id to connect paired probes (same as id/name for single probes).
--- 
--- Still have issues will true probe seq storage?  Some probes have mismatches and map multiple times therefore cannot use genome seq?
--- Seq storage not necessary for core, but maybe for FG analysis DB?
--- Have seq field/table?  Huge amount of probes and theoretically some could be >>1000bps.  If probes never overlapped would actually be less seq data as probe should be non-redundant within array.  Not the case, but may be easiest way to maintain true probe seq data.  Would negate the need for seq_region tables, this may hinder remapping. Then just have mismatches and cigar_line in probe_feature table instead of in both.  Would then also need build and chromosome tables.  How important is this?

--- As nimblegen probes are premapped to one location, can we capture this in probe_feature.analysis_id (e.g. vendor_mapping)?  Then alter this when we remap.  How can we maintain true match/mismatch pairs on remapping?  Only map PM probe, and then place MM probe at same location.  So mismatch in probe_feature would denote PM or MM probe only if probe/pair_index is valid and, vice versa.


--- Xref issue:  How are we going to consolidate vendor defined xrefs vs. ensembl core xrefs (e.g. affy) to ensure updating of xref table in core DB?


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
   `name` varchar(20) default NULL,
   `egroup_id` smallint(6) unsigned default NULL,
   `date` date NOT NULL default '0000-00-00',
   `primary_design_type` varchar(30) default NULL, 
   `description`  varchar(255) default NULL,
   PRIMARY KEY  (`experiment_id`),
   KEY `egroup_idx` (`egroup_id`),
   KEY `design_idx` (`primary_design_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- joint key on egroup and array?
--- design_type  = CHIP2 etc... (is also design type in ontology i.e. binding_site_identification)
--- Secondary design type may be redundant, so have associated_design_types table, containing MGED ontology types?  Or have ontology_types table to control input and have linker table with just IDs?  Too normalised?
--- epi_feature would also be redundant for non-CHIP2 experiments...have associated_target table?  Or just have target field, with "Experssion" for standard Xn chips?


--- Other fields:
---     egroup_member? (overkill?)

--- removed:
---    `array_id` int(11) unsigned default NULL, Cannot guarantee one array/chip set, especially as we may have to capture Nimblegen chips as individual chip sets.


--
-- Table structure for table `design_type`
--

DROP TABLE IF EXISTS `design_type`;
CREATE TABLE `design_type` (
   `design_type_id` int(11) unsigned NOT NULL auto_increment,
   `experiment_id` int(11) unsigned default NULL,	
   `name` varchar(40) default NULL,
   PRIMARY KEY  (`design_type_id`),
   KEY `channel_idx` (`experiment_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Handles design_types other than the primary experimental.primary_design_type
-- Do we need to duplicate primary_design_type in here?



--
-- Table structure for table `target`
--

DROP TABLE IF EXISTS `target`;
CREATE TABLE `target` (
   `target_id` int(11) unsigned NOT NULL auto_increment,
   `target_name` varchar(40) default NULL,
   `experiment_id` int(11) unsigned default NULL,
   `description`  varchar(255) default NULL,
   PRIMARY KEY  (`target_id`),
   KEY `experiment_idx` (`experiment_id`),	
   KEY `target_name_idx` (`target_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Indexes?
--- Change to associated_target table?
--- Table to contain Brno nomenclature (modification ontology?) 


--
-- Table structure for table `results`
--

DROP TABLE IF EXISTS `result`;
CREATE TABLE `result` (
   `result_id` int(11) unsigned NOT NULL auto_increment,
   `probe_id` int(11) unsigned default NULL,
   `score` double default NULL,
   `metric_id` int(11) unsigned default NULL,
   `channel_id` int(11) unsigned default NULL,
   PRIMARY KEY  (`result_id`),
   KEY `probe_idx` (`probe_id`),
   KEY `channel_idx` (`channel_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


--- REMOVEd experimental_id?   `experimental_chip_id` int(11) unsigned default NULL,
--- joint primary key with probe_feature_id?
--- metric default would be id for "RAW"..no, need to test has been specifically set, so NULL





--
-- Table structure for table `metric`
--

DROP TABLE IF EXISTS `metric`;
CREATE TABLE `metric` (
   `metric_id` int(11) unsigned NOT NULL auto_increment,
   `metric_name` varchar(20) default NULL,
   `description`  varchar(255) default NULL,
   PRIMARY KEY  (`metric_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--- Allows storage of none raw values
---Also needs to accommodate different normalisations 
--
-- Table structure for table `chip`
--

DROP TABLE IF EXISTS `experimental_chip`;
CREATE TABLE `experimental_chip` (
   `experimental_chip_id` int(11) unsigned NOT NULL auto_increment,
   `chip_unique_id` varchar(20) NOT NULL default '0',
   `experiment_id` int(11) unsigned default NULL,
   `array_chip_id` int(11) unsigned default NULL,
   PRIMARY KEY  (`experimental_chip_id`),
   KEY `experiment_idx` (`experiment_id`),
   KEY `chip_idx` (`chip_unique_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

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
   PRIMARY KEY  (`channel_id`),
   KEY `experimental_chip_idx` (`experimental_chip_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- type should be restricted to EXPERIMENTAL & CONTROL?
-- all other variables should be in experimental_variable table



--
-- Table structure for table `experimental_variable`
--

DROP TABLE IF EXISTS `experimental_variable`;
CREATE TABLE `experimental_variable` (
   `experimental_variable_id` int(11) unsigned NOT NULL auto_increment,
   `context_id` varchar(20) default NULL,	
   `name` varchar(40) default NULL,
   `unit` varchar(40) default NULL,
   `value` varchar(40) default NULL,
   PRIMARY KEY  (`experimental_variable_id`),
   KEY `channel_idx` (`context_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- context_id = table_name:table_name_id
-- Needs to accommodate QC flag and description, and reusing chips.
-- Have list of mandatory defs for each group, then also test MGED ontology and warn if not present.
-- status for chip?
-- Some of these are chip rather than channels specific
-- `tissue` varchar(40) NOT NULL default '',
-- `cell_line` varchar(40) NOT NULL default '',
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
-- Table structure for table `analysis`
--

DROP TABLE IF EXISTS `analysis`;
CREATE TABLE `analysis` (
   `analysis_id` int(11) unsigned NOT NULL auto_increment,
   `analysis_name` varchar(20) default NULL,
   `description` varchar(120) default NULL,
   PRIMARY KEY  (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


---Other fields as core.analysis?




--- Further thoughts:

--- More tables required for probe remapping? meta, rule tables etc?


--- Denormalise ---
--- To mart style to optimize queries? Mart style interface to export(to R)?
--- probe_set, probe, probe_feature and results could all e split on array vendor, format or experiment.design_type
---     probe > "array.vendor"_probe e.g. affy_probe, or "array.name"_probe e.g u133_probe, or "array.class"_probe e.g. CHIP2_probe


--- Mostly empty/unused fields which could be extracted to separate tables to reduce size?
---     cigar_line
