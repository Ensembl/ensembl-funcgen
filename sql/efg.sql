#
# Ensembl funcgen table definitions
# 


# Conventions:
#  - use lower case and underscores
#  - internal IDs are integers named tablename_id
#  - same name is given in foreign key relations

# DO NOT use '---' as a comment, this breaks in macosx mysql?!
# This also generally applies to all  --[^ ] !?

# Table documentation uses format parsed and defined here 


/** 
@header  Main feature tables
@desc    These define the various genomics features and their relevant associated tables.
*/


/**
@table  regulatory_feature
@desc   The table contains the features resulting from the regulatory build process.

@column regulatory_feature_id	Internal ID 
@column seq_region_id       	@link seq_region table ID
@column seq_region_start    	Start position of this feature
@column seq_region_end      	End position of this feature
@column seq_region_strand   	Strand orientation of this feature
@column display_label       	Text display label 
@column feature_type_id     	@link feature_type table ID
@column feature_set_id      	@link feature_set table ID
@column stable_id		      	Integer stable ID without ENSR prefix
@column bound_seq_region_start  Bound start position of this feature
@column bound_seq_region_end    Bound end position of this feature
@column binary_string			Binary representation for the underlying feature sets/types 
@column projected				Boolean, defines whether reg feat structure has been projected to this cell type 

@see seq_region
@see feature_set
@see feature_type
*/

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
  UNIQUE KEY `fset_seq_region_idx` (`feature_set_id`, `seq_region_id`,`seq_region_start`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `stable_id_idx` (`stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;



/**
@table  regulatory_attribute
@desc   Denormalised table defining links between a regulatory_feature and it's constituent 'attribute' features.

@column regulatory_feature_id	Internal ID 
@column attribute_feature_id	Table ID of attribute feature 
@column attribute_feature_table	Table name of attribute feature

@see annotated_feature
@see motif_feature
*/

DROP TABLE IF EXISTS `regulatory_attribute`;
CREATE TABLE `regulatory_attribute` (
  `regulatory_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated', 'external', 'motif') default NULL,
  PRIMARY KEY  (`regulatory_feature_id`, `attribute_feature_table`, `attribute_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;


-- remove external from enum

/**
@table  annotated_feature
@desc   Represents a genomic feature as the result of an analysis i.e. a ChIP or DNase1 peak call.

@column annotated_feature_id	Internal ID 
@column seq_region_id   	   	@link seq_region table ID
@column seq_region_start   		Start position of this feature
@column seq_region_end      	End position of this feature
@column seq_region_strand   	Strand orientation of this feature
@column display_label       	Text display label 
@column score               	Score derived from software
@column feature_set_id 	 		@link feature_set table ID
@column summit					Represents peak summit for those analyses which provide it (e.g. Swembl)

@see seq_region
@see feature_set
*/

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




/**
@table  motif_feature
@desc   The table contains genomic alignments of binding_matrix PWMs.

@column motif_feature_id    Primary key, internal ID 
@column binding_matrix_id   Foreign key to @link binding_matrix table
@column seq_region_id       Foreign key to @link seq_region table
@column seq_region_start    Start position of this feature
@column seq_region_end      End position of this feature
@column seq_region_strand   Strand orientation of this feature
@column display_label       Text display label 
@column score               Score derived from alignment software (e.g.MOODS)
@column interdb_stable_id   Unique key, provides linkability between DBs       

@see associated_motif_feature
@see regulatory_attribute
*/


DROP TABLE IF EXISTS `motif_feature`;
CREATE TABLE `motif_feature` (
  `motif_feature_id` int(10) unsigned NOT NULL auto_increment,
  `binding_matrix_id` INT(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) default NULL,
  `score` double default NULL,	
  `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL,
  PRIMARY KEY  (`motif_feature_id`),
  UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `binding_matrix_idx` (`binding_matrix_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  associated_motif_feature
@desc   The table provides links between motif_features and annotated_features representing peaks of the relevant transcription factor.

@column annotated_feature_id    @link annotated_feature table ID
@column motif_feature_id    	@link motif_feature table ID

@see associated_motif_feature
@see regulatory_attribute
*/


DROP TABLE IF EXISTS `associated_motif_feature`;
CREATE TABLE `associated_motif_feature` (
   `annotated_feature_id` int(10) unsigned NOT NULL,
   `motif_feature_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`annotated_feature_id`, `motif_feature_id`),
   KEY `motif_feature_idx` (`motif_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  binding_matrix
@desc   Contains information defining a specific binding matrix(PWM) as defined by the linked analysis e.g. Jaspar.

@column binding_matrix_id			Internal table ID
@column name			Name of PWM
@column feature_type_id			@link feature_type table ID.
@column frequencies			Matrix defining frequencing for each base at each position
@column description			Text description
@column analysis_id			@link analysis table ID
@column threshold			Minimum score for Motif Features for this matrix

@see motif_feature
@see feature_type
@see analysis
*/

DROP TABLE IF EXISTS `binding_matrix`;
CREATE  TABLE `binding_matrix` (
 `binding_matrix_id` INT(10) unsigned NOT NULL auto_increment,
 `name` VARCHAR(45) NOT NULL,
 `feature_type_id` int(10) unsigned NOT NULL,
 `frequencies` VARCHAR(1000) NOT NULL,
 `description` VARCHAR(255) NULL,
 `analysis_id` smallint(5) unsigned NOT NULL,
 `threshold` double default NULL,
 PRIMARY KEY (`binding_matrix_id`) ,
 KEY `feature_type_idx` (`feature_type_id`),
 KEY `name_analysis_idx` (`name`, `analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




/**
@table  external_feature
@desc   The table contains imports from externally curated resources e.g. cisRED, miRanda, VISTA, redFLY etc.

@column external_feature_id		Internal ID 
@column seq_region_id			@link seq_region table ID
@column seq_region_start 		Start position of this feature
@column seq_region_end      	End position of this feature
@column seq_region_strand   	Strand orientation of this feature
@column display_label       	Text display label 
@column feature_type_id     	@link feature_type table ID
@column feature_set_id      	@link feature_set table ID

@see seq_region
@see feature_set
@see feature_type
*/

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
  `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL,	
  PRIMARY KEY  (`external_feature_id`),
  UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


/** 
@table  result_feature
@desc   Represents the mapping of a raw/normalised signal. This is optimised for the web display
		in two ways:<br>
	&nbsp;&nbsp;&nbsp;&nbsp;1 Data compression by collection into different sized windows or bins.<br><br>
	&nbsp;&nbsp;&nbsp;&nbsp;2 For array data it also provides an optimised view of a probe_feature and associated result.
			
@column result_feature_id	Internal ID
@column result_set_id		@link result_set table ID 
@column seq_region_id       Foreign key to @link seq_region table
@column seq_region_start    Start position of this feature
@column seq_region_end      End position of this feature
@column seq_region_strand   Strand orientation of this feature
@column window_size			Size of window in base pairs
@column scores				BLOB of window scores for this region 

@see result_set
@see seq_region
*/

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

-- 0 window_size partition optimised for natural resolution of arrays
-- The rest are optimised for the default levels of zoom, plus a few intermediates 
-- No primary key as we need to parition on window_size, which would not be the whole of the primary key
-- set_window_region_idx can't be UNIQUE as there may be duplicates in the 0 wsize collections
-- i.e. two or more probe features with the same start, originating from replicate probes
-- or a probe seq which is a substr of another probe.
-- Remove partitions and limit to 0 window_size data, as others are now in DBFiles



/**
@table  probe_feature
@desc   The table contains genomic alignments @link probe entries.

@column probe_feature_id	Internal ID 
@column seq_region_id       Foreign key to @link seq_region table
@column seq_region_start    Start position of this feature
@column seq_region_end      End position of this feature
@column seq_region_strand   Strand orientation of this feature
@column probe_id			@link probe table ID
@column analysis_id			@link analysis table ID
@column mismatches			Integer, the number of bp matches for this alignment
@column cigar_line			Extended cigar line format representation of the alignment as defined here http://samtools.sourceforge.net/SAM-1.3.pdf.<br>In summary:
<br>&nbsp;&nbsp;&nbsp;&nbsp;= Seq/Alignment Match
<br>&nbsp;&nbsp;&nbsp;&nbsp;M Alignment match/Seq mismatch
<br>&nbsp;&nbsp;&nbsp;&nbsp;X Seq/Alignment mismatch
<br>&nbsp;&nbsp;&nbsp;&nbsp;D Deletion
<br>&nbsp;&nbsp;&nbsp;&nbsp;S Soft clipping, used for overhanging cdna alignments where genomic seq is unknown  

@see probe
@see seq_region
@see analysis
*/

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



/** 

@header  Set tables
@desc    Sets are containers for distinct sets of raw and/or processed data.

*/

/** 
@table  data_set
@desc   Defines highest level data container for associating the result of an analysis and the input data to that analysis e.g. Seq alignments(Input/ResultSet) and peak calls (FeatureSet)

@column data_set_id		Internal ID
@column feature_set_id	Product @link feature_set table ID
@column name            Name of data set

@see feature_set
@see supporting_set
*/

DROP TABLE IF EXISTS `data_set`;
CREATE TABLE `data_set` (
   `data_set_id` int(10) unsigned NOT NULL auto_increment,
   `feature_set_id` int(10) unsigned default '0',
   `name` varchar(100) default NULL,
   PRIMARY KEY  (`data_set_id`, `feature_set_id`),
   UNIQUE KEY `name_idx` (name)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




/** 
@table  supporting_set
@desc   Defines association between @link data_set and underlying/supporting data.

@column data_set_id			Internal ID
@column supporting_set_id	Table ID of supporting set
@column type				Type of supporting set e.g. result, feature or input set.

@see data_set
@see result_set
@see feature_set
@see input_set
*/


DROP TABLE IF EXISTS `supporting_set`;
CREATE TABLE `supporting_set` (
   `data_set_id` int(10) unsigned NOT NULL,
   `supporting_set_id` int(10) unsigned NOT NULL,
   `type` enum('result','feature','input') default NULL,
   PRIMARY KEY  (`data_set_id`, `supporting_set_id`),
   KEY `type_idx` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/** 
@table  feature_set
@desc   Container for genomic features defined by the result of an analysis e.g. peaks calls or regulatory features.

@column feature_set_id		Internal ID
@column feature_type_id		Table ID for @link feature_type
@column analysis_id			Table ID for @link analysis
@column cell_type_id		Table ID for @link cell_type
@column name				Name for this feature set
@column type				Type of features contained e.g. annotated, external or regualtory
@column description			Text description
@column display_label		Shorter more readable version of name

@see data_set
@see cell_type
@see feature_type
@see analysis
@see annotated_feature
@see regulatory_feature
@see external_feature
*/

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
   `experiment_id` int(10) unsigned default '0',
   PRIMARY KEY  (`feature_set_id`),
   KEY `feature_type_idx` (`feature_type_id`),
   KEY `experiment_idx` (`experiment_id`),
   UNIQUE KEY `name_idx` (name)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/** 
@table  result_set
@desc   Container for raw/signal data, used as input to an analysis or for visualisation of the raw signal i.e. a wiggle track.

@column result_set_id		Internal ID
@column analysis_id			Table ID for @link analysis
@column name				Name for this feature set
@column feature_type_id		Table ID for @link feature_type
@column cell_type_id		Table ID for @link cell_type

@see result_set_input 
@see cell_type
@see feature_type
@see analysis
@see supporting_set
@see dbfile_registry
*/

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


/** 
@table  result_set_input
@desc   Link table between @link result_set and it's contstituents which can vary between an array experiment (experimental_chip / channel) and a sequencing experiment (input_set). Note the joint primary key as inputs can be re-used between result sets.

@column result_set_input_id		Internal ID
@column result_set_id			@link result_set table ID 
@column table_id				Table ID for input 
@column table_name				Table name for input e.g. @link input_set, @link experimental_chip, @link channel

@see result_set
@see input_set
@see experimental_chip
@see channel
*/

DROP TABLE IF EXISTS `result_set_input`;
CREATE TABLE `result_set_input` (
   `result_set_input_id` int(10) unsigned NOT NULL auto_increment,
   `result_set_id` int(10) unsigned NOT NULL,
   `table_id` int(10) unsigned NOT NULL,
   `table_name` enum('experimental_chip','channel','input_set') DEFAULT NULL,
   PRIMARY KEY  (`result_set_input_id`, `result_set_id`),
   UNIQUE KEY `rset_table_idname_idx` (`result_set_id`, `table_id`, `table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



/**
@table  dbfile_registry

@desc   This generic table contains a simple registry of paths to support 
        flat file (DBFile) access. This should be left joined from the
        relevant adaptor e.g. ResultSetAdaptor.

@column table_id   Primary key of linked dbfile entity e.g. @link result_set or @link analysis
@column table_name Name of linked table
@column path       Either a full filepath or a directory which the API will use to build the filepath

@see result_set
*/


DROP TABLE IF EXISTS `dbfile_registry`;
CREATE TABLE `dbfile_registry` (
   `table_id` int(10) unsigned NOT NULL,
   `table_name` varchar(32)NOT NULL,	
   `path` varchar(255) NOT NULL,
   PRIMARY KEY  (`table_id`, `table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



/** 
@table  input_set
@desc   Defines a distinct set input data which is not imported into the DB, but used for some analysis e.g. a BAM file.

@column input_set_id		Internal ID
@column experiment_id		Table ID for @link experiment
@column feature_type_id		Table ID for @link feature_type
@column cell_type_id		Table ID for @link cell_type
@column format				Format of the experiment e.g. SEQUENCING
@column vendor				Name of technology vendor if appropriate
@column name				Name of input_set
@column type				Type of feature imported as result of analysis e.g. result, annotated

@see experiment
@see cell_type
@see feature_type
@see result_set_input
@see supporting_set
*/

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



/** 
@table  input_subset
@desc   Defines a file from an input_set, required for import tracking and recovery.

@column input_subset_id		Internal ID
@column input_set_id		@link input_set table ID
@column name				Name of input_subset e.g. file path

@see input_set
*/

DROP TABLE IF EXISTS `input_subset`;
CREATE TABLE `input_subset` (
   `input_subset_id` int(10) unsigned NOT NULL auto_increment,
   `input_set_id` int(10) unsigned NOT NULL,
   `name` varchar(100) NOT NULL, -- filename?	
   PRIMARY KEY  (`input_subset_id`), 
   UNIQUE KEY `set_name_dx` (`input_set_id`, `name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=30;



/** 
@header Array design tables
*/


/**
@table  array
@desc   Contains information defining an array or array set.

@column array_id	Internal ID
@column name     	Name of array
@column format      Format of array e.g. EXPRESSION, TILED
@column vendor		Name of array vendor e.g. AFFY
@column description Text description
@column type		Array type e.g. OLIGO, PCR
@column class		Array class e.g. AFFY_ST, ILLUMINA_INFINIUM  

@see    array_chip
*/

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

-- Could clean up/rename format, type and class


/**
@table  array_chip
@desc   Represents the individual array chip design as part of an array or array set.

@column array_chip_id	Internal ID
@column design_id  	    ID/Accession defined by vendor
@column array_id		Table ID from @link array
@column name     		Name of array_chip

@see array
@see probe
@see experimental_chip
*/

DROP TABLE IF EXISTS `array_chip`;
CREATE TABLE `array_chip` (
   `array_chip_id` int(10) unsigned NOT NULL auto_increment,
   `design_id` varchar(20) default NULL,
   `array_id` int(10) unsigned NOT NULL,
   `name` varchar(40) default NULL,
    PRIMARY KEY  (`array_chip_id`),
   UNIQUE KEY `array_design_idx` (`array_id`, `design_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




/**
@table  probe_set
@desc   The table contains information about probe sets.

@column probe_set_id		Internal ID 
@column name				Name of the probe set
@column size				Integer size of the probe set i.e. how many probe is contains
@column family				Generic descriptor for probe_set e.g. ENCODE_REGIONS, RANDOM etc. Currently not used

@see probe
*/

DROP TABLE IF EXISTS `probe_set`;
CREATE TABLE `probe_set` (
   `probe_set_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(100) NOT NULL,
   `size` smallint(6) unsigned NOT NULL,
   `family` varchar(20) default NULL,
   PRIMARY KEY  (`probe_set_id`),
	KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Clashing names may(unlikely) exist between array vendors, this is currently allowed
-- and handled with the API via the probe associations


/**
@table  probe
@desc   Defines individual probe designs across one or more array_chips. Note: The probe sequence is not stored. 

@column probe_id			Tnternal ID 
@column probe_set_id		@link probe_set table_id
@column name				Name of the probe set
@column length				Integer bp length of the probe
@column array_chip_id		@link array_chip table_id
@column class				Class of the probe e.g. CONTROL, EXPERIMENTAL etc.
@column description			Text description

@see probe_feature
@see probe_set
@see array_chip 
*/

DROP TABLE IF EXISTS `probe`;
CREATE TABLE `probe` (
   `probe_id` int(10) unsigned NOT NULL auto_increment,
   `probe_set_id` int(10) unsigned default NULL,
   `name` varchar(100) NOT NULL,
   `length` smallint(6) unsigned NOT NULL,
   `array_chip_id` int(10) unsigned NOT NULL,
   `class` varchar(20) default NULL,
   `description` varchar(255) DEFAULT NULL,
    PRIMARY KEY  (`probe_id`, `name`, `array_chip_id`),
    KEY `probe_set_idx`  (`probe_set_id`),
    KEY `array_chip_idx` (`array_chip_id`),
	KEY `name_idx`       (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  probe_design
@desc   Stores data from array design analyses.

@column probe_id			Internal ID 
@column analysis_id			@link analysis table ID
@column coord_system_id		@link coord_system table ID
@column score				Double analysis score value

@see probe
@see analysis
@see coord_system
*/

DROP TABLE IF EXISTS `probe_design`;
CREATE TABLE `probe_design` (
   `probe_id` int(10) unsigned NOT NULL,
   `analysis_id` smallint(5) unsigned NOT NULL,
   `score` double default NULL,	
   `coord_system_id` int(10) unsigned NOT NULL,
    PRIMARY KEY  (`probe_id`, `analysis_id`, `coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;




/** 
@header	Experiment tables
@desc	These define the experimental meta and raw data .
*/


/**
@table  experiment
@desc   Stores data high level meta data about individual experiments

@column experiment_id			Internal ID 
@column name			Name of experiment
@column experimental_group_id			@link experimental_group table ID
@column date			Date of experiment
@column primary_design_type			e.g. binding_site_identification, preferably EFO term
@column archive_id			ENA experiment identifier enabling access to specific raw data
@column data_url			When no accession_id exists, an alternative url to get the data 
@column description			Text description
@column	mage_xml_id			@link mage_xml table_id for array experiments

@see experimental_group
@see mage_xml
@see experimental_chip
@see input_set 
*/

DROP TABLE IF EXISTS `experiment`;
CREATE TABLE `experiment` (
   `experiment_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(100) default NULL,
   `experimental_group_id` smallint(6) unsigned default NULL,
   `date` date default '0000-00-00',
   `primary_design_type` varchar(30) default NULL, 
   `description`  varchar(255) default NULL,
   `mage_xml_id` int(10) unsigned default NULL,
   `archive_id` varchar(20) default NULL, 
   `data_url` varchar(255) default NULL, 
   PRIMARY KEY  (`experiment_id`),
   UNIQUE KEY `name_idx` (`name`),
   UNIQUE KEY `archive_idx` (`archive_id`),
   KEY `design_idx` (`primary_design_type`),
   KEY `experimental_group_idx` (`experimental_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Can probably remove now date (and primary_design_type?) as we don't support this level of meta data
-- also have associated_design_types/experimental_design table, containing MGED/EFO ontology types? 
-- change primary_design_type to design_type_id


/**
@table  experimental_group

@desc   Contains experimental group info i.e. who produced data sets.

@column experimental_group_id  Internal ID
@column name                   Name of group
@column location               Geographic location of group
@column contact                Contact details e.g. email
@column url                    Url for Project page
@column description            Text description

@see experiment
*/

DROP TABLE IF EXISTS `experimental_group`;
CREATE TABLE `experimental_group` (
   `experimental_group_id` smallint(6) unsigned NOT NULL auto_increment,
   `name` varchar(40) NOT NULL,
   `location` varchar(120) default NULL,
   `contact` varchar(40) default NULL,
   `description` varchar(255) default NULL,
   `url` varchar(255) default NULL,	
   `is_project` boolean default False,
   PRIMARY KEY  (`experimental_group_id`),
   UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



/**
@table  mage_xml
@desc   Contains MAGE-XML for array based experiments.

@column mage_xml_id		Internal table ID
@column xml				XML text field

@see experiment
*/

DROP TABLE IF EXISTS `mage_xml`;
CREATE TABLE `mage_xml` (
   `mage_xml_id` int(10) unsigned NOT NULL auto_increment,
   `xml` text,
   PRIMARY KEY  (`mage_xml_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  experimental_chip
@desc   Represents the physical instance of an @link array_chip used in an @link experiment.

@column experimental_chip_id	Internal ID
@column unique_id 		 	    Unique ID assigned by vendor
@column experiment_id			@link experiment table ID
@column array_chip_id     		@link array_chip table ID
@column feature_type_id    		@link feature_type table ID
@column cell_type_id     		@link cell_type table ID
@column biological_replicate	Name of biological replicate
@column technical_replicate		Name of technical replicate

@see array_chip 
@see experiment 
@see result_set_input 
@see channel 
@see feature_type
@see cell_type
*/

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

-- Can't implement unique key for unique_id as there may be clashes between vendors
-- and chips could potentially be re-used in another exp


/**
@table  channel
@desc   Represents an individual channel from an experimental_chip.

@column channel_id				Internal ID
@column experimental_chip_id	@link external_chip table ID
@column sample_id 		 	    Sample ID
@column dye						Name of dye used for this channel e.g. Cy3, Cy5
@column type					Type of channel i.e. EXPERIMENTAL or TOTAL (input)

@see  experimental_chip
*/

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


/** 
@table  result
@desc   Contains a score or intensity value for an associated probe location on a particular experimental_chip.

@column result_id				Internal ID
@column probe_id				@link probe table ID 
@column score					Intensity value (raw or normalised)
@column result_set_input_id		@link result_set_input table ID 
@column X						X coordinate of probe location on experimental_chip
@column Y						Y coordinate of probe location on experimental_chip

@see probe
@see result_set_input
*/

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



/** 

@header  Ancilliary tables
@desc    These contain data types which are used across many of the above tables and are quite often denormalised to store generic associations to several table, this avoids the need for multiple sets of similar tables.

*/

/** 
@table  feature_type

@desc   Contains information about different types/classes of feature e.g. Brno nomenclature, Transcription Factor names etc.

@column feature_type_id   Primary key, internal ID
@column name              Name of feature_type
@column class             Class of feature_type
@column description       Text description
@column so_accession      Sequence ontology accession
@column so_name           Sequence ontology name

@see associated_feature_type 
@see feature_set 
@see result_set 
@see input_set 
@see experimental_chip 
@see regulatory_feature 
@see external_feature
*/

DROP TABLE IF EXISTS `feature_type`;
CREATE TABLE `feature_type` (
	`feature_type_id` int(10) unsigned NOT NULL auto_increment,
	`name` varchar(40) NOT NULL,
	`class` enum('Insulator', 'DNA', 'Regulatory Feature', 'Histone', 'RNA', 
				'Polymerase', 'Transcription Factor', 'Transcription Factor Complex', 
				'Regulatory Motif',  'Enhancer', 'Expression', 'Pseudo', 
				'Open Chromatin', 'Search Region', 'Association Locus') default NULL,
	`description`  varchar(255) default NULL,
	`so_accession` varchar(64) DEFAULT NULL,
 	`so_name` varchar(255) DEFAULT NULL,
   PRIMARY KEY  (`feature_type_id`),
   UNIQUE KEY `name_class_idx` (`name`, `class`),
   KEY `so_accession_idx` (`so_accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  associated_feature_type
@desc   Link table providing many to many mapping for feature_type entries.

@column table_id			Internal table_id of linked table
@column table_name     	 	Name of linked table
@column feature_type_id		Internal table_id of linked @link feature_type

@see  feature_type
*/



DROP TABLE IF EXISTS `associated_feature_type`;
CREATE TABLE `associated_feature_type` (
   `table_id` int(10) unsigned NOT NULL,
   `table_name` enum('annotated_feature', 'external_feature', 'regulatory_feature', 'feature_type') NOT NULL,
   `feature_type_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`table_id`, `table_name`, `feature_type_id`),
   KEY `feature_type_index` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



/** 
@table  cell_type

@desc   Contains information about cell/tissue types.

@column cell_type_id	Internal ID
@column name        	Name of cell/tissue
@column display_label	Short display label
@column description 	Text description
@column gender			Gender i.e. male or female
@column efo_id          Experimental Factor Ontology ID

@see feature_set
@see result_set
@see input_set
@see experimental_chip
*/

DROP TABLE IF EXISTS `cell_type`;
CREATE TABLE `cell_type` (
   `cell_type_id` int(10) unsigned NOT NULL auto_increment,
   `name`  varchar(120) not NULL,
   `display_label` varchar(20) default NULL,
   `description` varchar(80) default NULL,
   `gender` enum('male', 'female') default NULL,
   `efo_id` varchar(20) DEFAULT NULL,
   PRIMARY KEY  (`cell_type_id`),
   UNIQUE KEY `name_idx` (`name`),
   UNIQUE KEY `efo_idx` (`efo_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- bto?

/**
@table  experimental_design
@desc   Denormalised link table to allow many to many design_type associations.

@column design_type_id			@link design_type table ID
@column table_name				Name of linked table
@column table_id				Table ID of linked record

@see design_type
*/

DROP TABLE IF EXISTS `experimental_design`;
CREATE TABLE `experimental_design` (
   `design_type_id` int(10) unsigned NOT NULL auto_increment,
   `table_name` varchar(40) default NULL,
   `table_id` int(10) unsigned default NULL,	
   PRIMARY KEY  (`design_type_id`, `table_name`, `table_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- rename associated_design_type, or replace with ontology xrefs?


/**
@table  design_type
@desc   Contains extra information about experimental designs, preferably ontology terms. 

@column design_type_id	Internal ID
@column name			Name of design type

@see experimental_design
*/

-- extend to support ontologies properly or remove

DROP TABLE IF EXISTS `design_type`;
CREATE TABLE `design_type` (
   `design_type_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(255) default NULL,
   PRIMARY KEY  (`design_type_id`),
   KEY `design_name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- add description?


/**
@table  status
@desc   Denormalised table associating funcgen records with a status.

@column table_id		Table ID of associated record
@column table_name		Table name of associated record
@column status_name_id	@link status_name table ID


@see status_name
*/

DROP TABLE IF EXISTS `status`;
CREATE TABLE `status` (
   `table_id` int(10) unsigned default NULL,
   `table_name` varchar(32) default NULL,	
   `status_name_id` int(10) unsigned NOT NULL,
   PRIMARY KEY  (`table_id`, `table_name`, `status_name_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- 32 is max length for mysql table name


/**
@table  status_name
@desc   Simple table to predefine name of status.

@column status_name_id	Internal ID
@column name			Name of status e.g. IMPORTED, DISPLAYBLE etc.

@see status
*/

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






/** 

@header  Core tables
@desc    These are exact clones of the corresponding core schema tables. See <a href='../core/core_schema.html'>core schema docs</a> for more details.

*/

/**
@table analysis
@desc Usually describes a program and some database that together are used to create a feature on a piece of sequence.
Each feature is marked with an analysis_id. The most important column is logic_name, which is used by the webteam to render a feature
 correctly on contigview (or even retrieve the right feature).
Logic_name is also used in the pipeline to identify the analysis which has to run in a given status of the pipeline.
The module column tells the pipeline which Perl module does the whole analysis, typically a RunnableDB module.

@column analysis_id                 Internal ID
@column created                     Date to distinguish newer and older versions off the same analysis.
@column logic_name                  String to identify the analysis. Used mainly inside pipeline.
@column db                          Database name.
@column db_version                  Database version.
@column db_file                     File system location of the database.
@column program                     The binary used to create a feature. 
@column program_version             The binary version.
@column program_file                File system location of the binary.
@column parameters                  A parameter string which is processed by the perl module.
@column module                      Perl module names (RunnableDBS usually) executing this analysis.
@column module_version              Perl module version.
@column gff_source                  How to make a gff dump from features with this analysis.
@column gff_feature                 How to make a gff dump from features with this analysis.

@see analysis_description

*/


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


/**
@table analysis_description
@desc Allows the storage of a textual description of the analysis, as well as a "display label", primarily for the EnsEMBL web site.

@column analysis_id            Foreign key references to the @link analysis table.
@column description            Textual description of the analysis.
@column display_label          Display label for the EnsEMBL web site.
@column displayable            Flag indicating if the analysis description is to be displayed on the EnsEMBL web site.
@column web_data               Other data used by the EnsEMBL web site.

@see analysis

*/


DROP TABLE IF EXISTS `analysis_description`;
CREATE TABLE `analysis_description` (
  `analysis_id` smallint(5) unsigned NOT NULL,
  `description` text,
  `display_label` varchar(255) NOT NULL,
  `displayable` BOOLEAN NOT NULL default '1',
  `web_data` text,	
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table meta
@desc Stores data about the data in the current schema. Unlike other tables, data in the meta table is stored as key-value pairs. These data include details about the database, RegulatoryBuild and patches.  The species_id field of the meta table is used in multi-species databases and makes it possible to have species-specific meta key-value pairs.  The species-specific meta key-value pairs needs to be repeated for each species_id.  Entries in the meta table that are not specific to any one species, such as the schema.version key and any other schema-related information must have their species_id field set to NULL
. The default species_id, and the only species_id value allowed in single-species databases, is 1.


@column meta_id                    Internal identifier.
@column species_id                 Indentifies the species for multi-species databases.
@column meta_key                   Name of the meta entry, e.g. "schema_version".
@column meta_value                 Corresponding value of the key, e.g. "61".

*/


DROP TABLE IF EXISTS `meta`;
CREATE TABLE `meta` (
  `meta_id` int(10) NOT NULL auto_increment,
  `species_id` int(10) unsigned default '1',
  `meta_key` varchar(46) NOT NULL,
  `meta_value` varchar(950) NOT NULL,
  PRIMARY KEY  (`meta_id`),
  UNIQUE KEY `species_key_value_idx` (`species_id`,`meta_key`,`meta_value`),
  KEY `species_value_idx` (`species_id`,`meta_value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- Add empty schema_version entry for patches to update
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, "schema_version", "");
INSERT INTO meta (meta_key, meta_value) VALUES ('schema_type', 'funcgen');

-- Also need to add the current release patch entries here?


/**
@table meta_coord
@desc Describes which co-ordinate systems the different feature tables use.

@column table_name              Ensembl database table name.
@column coord_system_id         Table ID for @link coord_system
@column max_length              Longest sequence length.    

@see coord_system

*/

DROP TABLE IF EXISTS `meta_coord`;
CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `max_length` int(11) default NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table identity_xref
@desc Describes how well a particular xref object matches the EnsEMBL object.

@column object_xref_id        Foreign key references to the @link object_xref table.
@column xref_identity         Percentage identity.
@column ensembl_identity      Percentage identity.
@column xref_start            Xref sequence start.
@column xref_end              Xref sequence end.
@column ensembl_start         Ensembl sequence start.
@column ensembl_end           Ensembl sequence end.
@column cigar_line            Used to encode gapped alignments.
@column score                 Match score.
@column evalue                Match evalue.

@see object_xref

*/

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




/**
@table external_synonym
@desc Some xref objects can be referred to by more than one name. This table relates names to xref IDs.

@column xref_id           Foreign key references @link xref table
@column synonym           Synonym

@see xref

*/


DROP TABLE IF EXISTS external_synonym;
CREATE TABLE external_synonym (
  xref_id                     INT(10) UNSIGNED NOT NULL,
  synonym                     VARCHAR(100) NOT NULL, 
  PRIMARY KEY (xref_id, synonym),
  KEY name_index (synonym)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=20;


/**
@table external_db
@desc Stores data about the external databases in which the objects described in the xref table are stored.

@column external_db_id              Internal identifier.
@column db_name                     Database name.     
@column db_release                  Database release.      
@column status                      Status, e.g. 'KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO'. 
@column dbprimary_acc_linkable      Indicates if primary a accession can be linked to from the EnsEMBL web site.
@column priority                    Determines which one of the xrefs will be used as the gene name.
@column db_display_name             Database display name.
@column type                        Type, e.g. 'ARRAY', 'ALT_TRANS', 'ALT_GENE', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL'.
@column secondary_db_name           Secondary database name.
@column secondary_db_table          Secondary database table.
@column description                 Description.


@see xref
@see unmapped_object

*/


DROP TABLE IF EXISTS external_db;
CREATE TABLE external_db (
  external_db_id 	          SMALLINT(5) UNSIGNED NOT NULL auto_increment,
  db_name                     VARCHAR(100) NOT NULL,
  db_release                  VARCHAR(255),
  status                      ENUM('KNOWNXREF','KNOWN','XREF','PRED','ORTH', 'PSEUDO') NOT NULL,
  dbprimary_acc_linkable      BOOLEAN DEFAULT 1 NOT NULL,
  priority                    INT NOT NULL,
  db_display_name             VARCHAR(255),
  type                        ENUM('ARRAY', 'ALT_TRANS', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL') default NULL,
  secondary_db_name           VARCHAR(255) DEFAULT NULL,
  secondary_db_table          VARCHAR(255) DEFAULT NULL,
  description                 TEXT,
  PRIMARY KEY (external_db_id) ,
  UNIQUE KEY db_name_release_idx (db_name, db_release)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=80;



/**
@table ontology_xref
@desc This table associates ontology terms/accessions to Ensembl objects (primarily EFO/SO). NOTE: Currently not in use

@column object_xref_id          Foreign key references to the @link object_xref table.
@column source_xref_id          Foreign key references to the @link xref table.
@column linkage_type            Defines type of linkage

@see object_xref

*/


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

/**
@table unmapped_reason
@desc Describes the reason why a mapping failed.

@column unmapped_reason_id           Internal identifier.
@column summary_description          Summarised description.
@column full_description             Full description.

@see unmapped_object
*/

DROP TABLE if EXISTS unmapped_reason;
CREATE TABLE `unmapped_reason` (
  `unmapped_reason_id` smallint(5) unsigned NOT NULL auto_increment,
  `summary_description` varchar(255) default NULL,
  `full_description` varchar(255) default NULL,
  PRIMARY KEY  (`unmapped_reason_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/** 

@header  Core like tables
@desc    These are almost exact clones of the corresponding core schema tables. Some contain extra fields or different enum values to support the funcgen schema

*/

/**
@table xref
@desc Holds data about objects which are external to EnsEMBL, but need to be associated with EnsEMBL objects.
Information about the database that the external object is stored in is held in the external_db table entry referred to by the external_db column.

@column xref_id                 Internal identifier.
@column external_db_id          Foreign key references to the @link external_db table.
@column dbprimary_acc           Primary accession number.
@column display_label           Display label for the EnsEMBL web site.
@column version                 Object version.
@column description             Object description.
@column info_type               Class of the xref information e.g. CODING
@column info_text               Text

@see 
@see external_db
@see external_synonym
@see xref

*/

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


/**
@table object_xref
@desc Describes links between Ensembl objects and objects held in external databases.  The Ensembl object can be one of several types; the type is held in the ensembl_object_type column.  The ID of the particular Ensembl gene, translation or whatever is given in the ensembl_id column.  The xref_id points to the entry in the xref table that holds data about the external object.  Each Ensembl object can be associated with zero or more xrefs. An xref object can be associated with one or more Ensembl objects.

@column object_xref_id            Internal identifier.
@column ensembl_id                Foreign key references to the ensembl_object_type table e.g. @link probe_set
@column ensembl_object_type       Ensembl object type e.g ProbeSet etc.
@column xref_id                   Foreign key references to the @link xref table.
@column linkage_annotation        Additional annotation on the linkage.
@column analysis_id               Foreign key references to the @link analysis table.

@see xref
@see identity_xref

*/


DROP TABLE IF EXISTS `object_xref`;
CREATE TABLE object_xref (
  object_xref_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  ensembl_id                  INT(10) UNSIGNED NOT NULL, 
  ensembl_object_type         ENUM('Experiment', 'RegulatoryFeature', 'ExternalFeature', 'AnnotatedFeature', 'FeatureType', 'ProbeSet', 'Probe', 'ProbeFeature') not NULL,
  xref_id                     INT UNSIGNED NOT NULL,
  linkage_annotation          VARCHAR(255) DEFAULT NULL,
  analysis_id                 SMALLINT(5) UNSIGNED NOT NULL,
  PRIMARY KEY (`object_xref_id`),
  UNIQUE KEY `xref_idx` (`xref_id`,`ensembl_object_type`,`ensembl_id`,`analysis_id`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `ensembl_idx` (`ensembl_object_type`,`ensembl_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- Note we use case correct versions of object name to allow easy adaptor generation


/**
@table unmapped_object
@desc Describes why a particular external entity was not mapped to an ensembl one.

@column unmapped_object_id         Internal identifier.
@column type                       UnmappedObject type e.g. probe2transcript
@column analysis_id                Foreign key references to the @link analysis table.
@column external_db_id             Foreign key references to the @link external_db table.
@column identifier                 External database identifier.
@column unmapped_reason_id         Foreign key references to the @link unmapped_reason table.
@column query_score                Actual mapping query score.
@column target_score               Target mapping query score.
@column ensembl_id                 Foreign key references the table_if of the Ensembl object table e.g. @link probe_set
@column ensembl_object_type        Ensembl object type e.g. ProbeSet
@column parent                     Foreign key references to the @link dependent_xref table, in case the unmapped object is dependent on a primary external reference which wasn't mapped to an ensembl one. Not currently used for efg.

@see unmapped_reason
@see external_db
@see analysis
*/



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
  UNIQUE KEY `unique_unmapped_obj_idx` (`identifier`,`ensembl_id`,`parent`,`unmapped_reason_id`,`ensembl_object_type`,`external_db_id`),
  KEY `anal_exdb_idx` (`analysis_id`,`external_db_id`),
  KEY `id_idx` (`identifier`(50)),
  KEY `ext_db_identifier_idx` (`external_db_id`,`identifier`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table coord_system
@desc Stores information about the available co-ordinate systems for the species identified through the species_id field.
For each species, there must be one co-ordinate system that has the attribute "top_level" and one that has the attribute "sequence_level". 
NOTE: This has been extended from the core implementation to support multiple assemblies by referencing multiple core DBs.

@column coord_system_id      Internal identifier
@column name                 Co-oridinate system name, e.g. 'chromosome', 'contig', 'scaffold' etc. 
@column version              Assembly
@column rank                 Co-oridinate system rank
@column attrib               Co-oridinate system attrib (e.g. "top_level", "sequence_level")     
@column schema_build         Indentifies the schema_build version for the source core DB       
@column core_coord_system_id Table ID of the coord_system in the source core DB
@column species_id           Indentifies the species for multi-species databases       
@column is_current           This flags which coord_system entries are current with respect to mart(/website)

@see seq_region
@see meta_coord

*/


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

-- Do we need attrib and rank for eFG?


/**
@table seq_region
@desc Stores information about sequence regions from various core DBs.

@column seq_region_id            Internal identifier.
@column name                     Sequence region name.  
@column coord_system_id          Foreign key references to the @link coord_system table.   
@column core_seq_region_id       Table ID of the seq_region in the source core DB
@column schema_build             Indentifies the schema_build version for the source core DB       

@see coord_system
*/



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
-- hence schema_build denormalisation

-- Name is only required to enable us to add new seq_regions to the correct seq_region_id
-- It will never be used to retrieve a slice as we do that via the core DB
-- Usage: Pull back seq_region_id based schema_build and core_seq_region_id



