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

/**
* Ensembl funcgen table definitions
*
*
*
* Conventions:
*  - use lower case and underscores
*  - internal IDs are integers named tablename_id
*  - same name is given in foreign key relations
*
* DO NOT use '---' as a comment, this breaks in macosx mysql?!
* This also generally applies to all  --[^ ] !?
*
* Table documentation uses format parsed and defined here

*/

/**
@header  Main feature tables
@desc    These define the various genomics features and their relevant associated tables.
@colour  #FFCC66
@legend  #FFCC66 Main feature tables
*/


/**
@table  regulatory_feature
@desc   The table contains the features resulting from the regulatory build process.
@colour  #FFCC66

@column regulatory_feature_id   Internal ID
@column feature_type_id         @link feature_type ID
@column seq_region_id           seq_region ID
@column seq_region_start        Start position of this featurefeature_set
@column seq_region_end          End position of this feature
@column seq_region_strand       Strand orientation of this feature
@column stable_id               Integer stable ID without ENSR prefix *mnuhn: Not true, they do have this prefix*
@column bound_start_length      Distance between start of the feature and start of the bound region. Bound regions are used for promoters only. They define the flanking regions. It is an area that is predicted t
@column bound_end_length        Distance between end of the bound region and end of this feature
@column epigenome_count         Integer, number of epigenomes in which this feature is active
@column regulatory_build_id     @link regulatory_build ID

@see feature_type

*/

DROP TABLE IF EXISTS `regulatory_feature`;
CREATE TABLE `regulatory_feature` (
  `regulatory_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `stable_id` varchar(18) DEFAULT NULL,
  `bound_start_length` mediumint(3) unsigned NOT NULL,
  `bound_end_length` mediumint(3) unsigned NOT NULL,
  `epigenome_count` smallint(6) DEFAULT NULL,
  `regulatory_build_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`regulatory_feature_id`),
  -- The name "uniqueness_constraint_idx" is used in the
  -- RegulatoryFeatureAdaptor to catch issues regarding its violation.
  -- Changing the name means having to update it in the
  -- RegulatoryFeatureAdaptor as well.
  --
  UNIQUE KEY `uniqueness_constraint_idx` (`feature_type_id`,`seq_region_id`,`seq_region_strand`,`seq_region_start`,`seq_region_end`,`stable_id`,`bound_start_length`,`bound_end_length`,`regulatory_build_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `stable_id_idx` (`stable_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=latin1;

/**
@table  regulatory_evidence
@desc   Links a regulatory feature and the epigenome (via regulatory
        activity) to the underlying structure of epigenetic marks that the
        regulatory feature has in this epigenome.
@colour  #FFCC66

@column regulatory_activity_id  @link regulatory_activity
@column attribute_feature_id    Table ID of attribute feature
@column attribute_feature_table Table name of attribute feature

@see peak
@see regulatory_activity
*/

DROP TABLE IF EXISTS `regulatory_evidence`;
CREATE TABLE `regulatory_evidence` (
  `regulatory_activity_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated','motif') NOT NULL DEFAULT 'annotated',
  PRIMARY KEY (`regulatory_activity_id`,`attribute_feature_table`,`attribute_feature_id`),
  KEY `attribute_feature_idx` (`attribute_feature_id`,`attribute_feature_table`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  regulatory_activity
@desc   For every regulatory feature and epigenome that was a part of the
        regulatory build, this table links the regulatory feature to the
        predicted regulatory activity in this epigenome.
@colour  #FFCC66

@column regulatory_activity_id  Internal ID
@column regulatory_feature_id   @link regulatory_feature
@column activity                The predicted activity of the regulatory feature in the epigenome.
@column epigenome_id            @link epigenome

@see epigenome
@see regulatory_feature

*/

DROP TABLE IF EXISTS `regulatory_activity`;
CREATE TABLE `regulatory_activity` (
  `regulatory_activity_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_feature_id` int(10) unsigned DEFAULT NULL,
  `activity` enum('INACTIVE','REPRESSED','POISED','ACTIVE','NA') DEFAULT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`regulatory_activity_id`),
  UNIQUE KEY `uniqueness_constraint_idx` (`epigenome_id`,`regulatory_feature_id`),
  KEY `regulatory_feature_idx` (`regulatory_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  regulatory_build
@desc   Metadata for the regulatory build
@colour  #FFCC66

@column regulatory_build_id            Internal ID
@column name                           Name of the regulatory build
@column version                        Version of the regulatory build
@column initial_release_date           Date of initial release
@column last_annotation_update         Date of last annotation update
@column feature_type_id                @link feature_type
@column analysis_id                    @link analysis
@column is_current                     Set to true, if this entry refers to the current regulatory build
@column sample_regulatory_feature_id   @link regulatory_feature

@see feature_type
@see regulatory_feature
@see analysis

*/

DROP TABLE IF EXISTS `regulatory_build`;
CREATE TABLE `regulatory_build` (
  `regulatory_build_id` int(4) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(45) NOT NULL,
  `version` varchar(50) DEFAULT NULL,
  `initial_release_date` varchar(50) DEFAULT NULL,
  `last_annotation_update` varchar(50) DEFAULT NULL,
  `feature_type_id` int(4) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `is_current` tinyint(1) NOT NULL DEFAULT '0',
  `sample_regulatory_feature_id` int(10) unsigned,
  PRIMARY KEY (`regulatory_build_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  regulatory_build_epigenome
@desc   Table that links a regulatory build to the epigenomes that were used in it.
@colour  #FFCC66

@column regulatory_build_epigenome_id Internal ID
@column regulatory_build_id           @link regulatory_build
@column epigenome_id                  @link epigenome

@see regulatory_build
@see epigenome

*/

DROP TABLE IF EXISTS `regulatory_build_epigenome`;
CREATE TABLE `regulatory_build_epigenome` (
  `regulatory_build_epigenome_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_build_id` int(10) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`regulatory_build_epigenome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  segmentation_feature
@desc   Represents a genomic segment feature as the result of an segmentation
        analysis i.e. Segway or ChromHmm
@colour  #FFCC66

@column segmentation_feature_id Internal ID
@column feature_set_id          @link feature_set ID
@column feature_type_id         @link feature_type ID
@column seq_region_id           seq_region ID
@column seq_region_start        Start position of this feature
@column seq_region_end          End position of this feature
@column seq_region_strand       Strand orientation of this feature
@column score                   Score derived from software
@column display_label           Text display label

@see feature_set
@see feature_type
*/

DROP TABLE IF EXISTS `segmentation_feature`;
CREATE TABLE `segmentation_feature` (
  `segmentation_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `feature_set_id` int(10) unsigned DEFAULT NULL,
  `score` double DEFAULT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  PRIMARY KEY (`segmentation_feature_id`),
  UNIQUE KEY `fset_seq_region_idx` (`feature_set_id`,`seq_region_id`,`seq_region_start`),
  KEY `feature_type_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

/**
@table  segmentation_file
@desc   Table to store metadata about a segmentation file
@colour  #FFCC66

@column segmentation_file_id Internal ID
@column regulatory_build_id  @link regulatory_build
@column name                 A descriptive name of what is in the file.
@column analysis_id          @link analysis
@column epigenome_id         @link epigenome

@see regulatory_build
@see epigenome
@see analysis

*/

DROP TABLE IF EXISTS `segmentation_file`;
CREATE TABLE `segmentation_file` (
  `segmentation_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_build_id` int(4) unsigned DEFAULT NULL,
  `name` varchar(100) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`segmentation_file_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `epigenome_idx` (`epigenome_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  peak
@desc   Represents a genomic feature as the result of an analysis i.e. a ChIP or DNase1 peak call.
@colour  #FFCC66

@column peak_id    Internal ID
@column peak_calling_id         @link peak_calling ID
@column seq_region_id           seq_region ID
@column seq_region_start        Start position of this feature
@column seq_region_end          End position of this feature
@column seq_region_strand       Strand orientation of this feature
@column score                   Score derived from software
@column summit                  Represents peak summit for those analyses which provide it (e.g. Swembl)

@see peak_calling
*/

DROP TABLE IF EXISTS `peak`;
CREATE TABLE `peak` (
  `peak_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `score` double DEFAULT NULL,
  `peak_calling_id` int(10) unsigned NOT NULL,
  `summit` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`peak_id`),
  UNIQUE KEY `seq_region_feature_set_idx` (`seq_region_id`,`seq_region_start`,`peak_calling_id`),
  KEY `feature_set_idx` (`peak_calling_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=39;

/**
@table  peak_calling
@desc   Represents a peak calling analysis.
@colour  #FFCC66

@column peak_calling_id         Internal ID
@column name                    Name of the peak calling
@column display_label           Name for displaying on the website
@column feature_type_id         @link feature_type ID
@column analysis_id             @link analysis ID
@column alignment_id            @link alignment ID
@column epigenome_id            @link epigenome ID
@column experiment_id           @link experiment ID

@see feature_type
@see analysis
@see alignment
@see epigenome
@see experiment
*/

DROP TABLE IF EXISTS `peak_calling`;
CREATE TABLE `peak_calling` (
  `peak_calling_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(300) NOT NULL,
  `display_label` varchar(300) NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `alignment_id` int(10) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `experiment_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`peak_calling_id`),
  UNIQUE KEY `peak_calling_id_idx` (`peak_calling_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  motif_feature
@desc   The table contains genomic alignments of binding_matrix PWMs.
@colour  #FFCC66

@column motif_feature_id    Primary key, internal ID
@column binding_matrix_id   @link binding_matrix table
@column seq_region_id       seq_region table
@column seq_region_start    Start position of this feature
@column seq_region_end      End position of this feature
@column seq_region_strand   Strand orientation of this feature
@column display_label       Text display label
@column score               Score derived from alignment software (e.g.MOODS)
@column interdb_stable_id   Unique key, provides linkability between DBs

@see associated_motif_feature
@see binding_matrix
*/


DROP TABLE IF EXISTS `motif_feature`;
CREATE TABLE `motif_feature` (
  `motif_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `binding_matrix_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  `score` double DEFAULT NULL,
  `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL,
  PRIMARY KEY (`motif_feature_id`),
  UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `binding_matrix_idx` (`binding_matrix_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  mirna_target_feature
@desc   The table contains imports from externally curated resources e.g. cisRED, miRanda, VISTA, redFLY etc.
@colour  #FFCC66

@column mirna_target_feature_id Internal ID
@column feature_set_id          @link feature_set ID
@column feature_type_id         @link feature_type ID
@column seq_region_id           seq_region ID
@column accession               Accession number given by data source
@column display_label           Text display label
@column evidence                Evidence level provided by data source
@column interdb_stable_id       Unique key, provides linkability between DBs
@column method                  Method used to identify miRNA target
@column seq_region_start        Start position of this feature
@column seq_region_end          End position of this feature
@column seq_region_strand       Strand orientation of this feature
@column supporting_information  Additional information which does not fit another category

@see feature_set
@see feature_type
*/

DROP TABLE IF EXISTS `mirna_target_feature`;
CREATE TABLE `mirna_target_feature` (
  `mirna_target_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `feature_set_id` int(10) unsigned NOT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `accession` varchar(60) DEFAULT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  `evidence` varchar(60) DEFAULT NULL,
  `interdb_stable_id` int(10) unsigned DEFAULT NULL,
  `method` varchar(60) DEFAULT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `supporting_information` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`mirna_target_feature_id`),
  UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;


/**
@table  associated_motif_feature
@desc   The table provides links between motif_features and annotated_features representing peaks of the relevant transcription factor.
@colour  #FFCC66

@column annotated_feature_id    @link annotated_feature table ID
@column motif_feature_id        @link motif_feature table ID

*/

DROP TABLE IF EXISTS `associated_motif_feature`;
CREATE TABLE `associated_motif_feature` (
  `annotated_feature_id` int(10) unsigned NOT NULL,
  `motif_feature_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`annotated_feature_id`,`motif_feature_id`),
  KEY `motif_feature_idx` (`motif_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  binding_matrix
@desc   Contains information defining a specific binding matrix(PWM) as defined by the linked analysis e.g. Jaspar.
@colour  #FFCC66

@column binding_matrix_id  Internal table ID
@column analysis_id        @link analysis table ID
@column feature_type_id    @link feature_type table ID.
@column name               Name of PWM
@column frequencies        Matrix defining frequencing for each base at each position
@column description        Text description
@column threshold          Minimum score for Motif Features for this matrix

@see analysis
@see feature_type
*/

DROP TABLE IF EXISTS `binding_matrix`;
CREATE TABLE `binding_matrix` (
  `binding_matrix_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(45) NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  `frequencies` varchar(1000) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `threshold` double DEFAULT NULL,
  PRIMARY KEY (`binding_matrix_id`),
  UNIQUE KEY `name_analysis_idx` (`name`,`analysis_id`),
  KEY `feature_type_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  external_feature
@desc   The table contains imports from externally curated resources e.g. cisRED, miRanda, VISTA, redFLY etc.
@colour  #FFCC66

@column external_feature_id  Internal ID
@column feature_set_id       @link feature_set ID
@column feature_type_id      @link feature_type ID
@column seq_region_id        seq_region ID
@column seq_region_start     Start position of this feature
@column seq_region_end       End position of this feature
@column seq_region_strand    Strand orientation of this feature
@column display_label        Text display label
@column interdb_stable_id    Unique key, provides linkability between DBs

@see feature_set
@see feature_type
*/

DROP TABLE IF EXISTS `external_feature`;
CREATE TABLE `external_feature` (
  `external_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `feature_set_id` int(10) unsigned NOT NULL,
  `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL,
  PRIMARY KEY (`external_feature_id`),
  UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

/**
@table  external_feature_file
@desc   Table to store metadata about a file with features
@colour  #FFCC66

@column external_feature_file_id  Internal ID
@column name                      A name descriptive of the data in the file
@column analysis_id               @link analysis
@column epigenome_id              @link epigenome
@column feature_type_id           @link feature_type

@see data_file

*/

DROP TABLE IF EXISTS `external_feature_file`;
CREATE TABLE `external_feature_file` (
  `external_feature_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`external_feature_file_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `epigenome_idx` (`epigenome_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  probe_feature
@desc   The table contains genomic alignments @link probe entries.
@colour  #FFCC66

@column probe_feature_id    Internal ID
@column analysis_id         @link analysis table ID
@column probe_id            @link probe table ID
@column seq_region_id       seq_region table ID
@column seq_region_start    Start position of this feature
@column seq_region_end      End position of this feature
@column seq_region_strand   Strand orientation of this feature
@column mismatches          Integer, the number of bp matches for this alignment
@column cigar_line          Extended cigar line format representation of the alignment as defined here http://samtools.sourceforge.net/SAM-1.3.pdf.<br>
                            In summary:
                            <ul>
                              <li>= Seq/Alignment Match</li>
                              <li>M Alignment match/Seq mismatch</li>
                              <li>X Seq/Alignment mismatch</li>
                              <li>D Deletion</li>
                              <li>S Soft clipping, used for overhanging cdna alignments where genomic seq is unknown</li>
                            </ul>
@column hit_id              Id of the sequence on which the hit was initially made. Typically this will be the stable id of the transcript or the name of the sequence region.
@column source              The source of the sequence on which the probe was found. If set, this can be 'genomic' or 'transcript'

@see analysis
@see probe
*/

DROP TABLE IF EXISTS `probe_feature`;
CREATE TABLE `probe_feature` (
  `probe_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) NOT NULL,
  `seq_region_end` int(10) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `probe_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `mismatches` tinyint(4) NOT NULL,
  `cigar_line` varchar(50) DEFAULT NULL,
  `hit_id` varchar(255) DEFAULT NULL,
  `source` enum('genomic','transcript') DEFAULT NULL,
  PRIMARY KEY (`probe_feature_id`),
  KEY `probe_idx` (`probe_id`),
  KEY `seq_region_probe_probe_feature_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`,`probe_id`,`probe_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  probe_feature_transcript
@desc   The table maps probe_features to transcripts.
@colour  #FFCC66

@column probe_feature_transcript_id    Internal ID
@column probe_feature_id               @link probe_feature table ID
@column stable_id                      Stable id of the transcript to which it has been mapped
@column description                    Transcript description

@see probe_feature
*/

DROP TABLE IF EXISTS `probe_feature_transcript`;
CREATE TABLE `probe_feature_transcript` (
  `probe_feature_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_feature_id` int(10) unsigned DEFAULT NULL,
  `stable_id` varchar(128) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_feature_transcript_id`),
  KEY `probe_feature_transcript_id_idx` (`probe_feature_transcript_id`),
  KEY `probe_feature_id_idx` (`probe_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  feature_type
@desc   Contains information about different types/classes of feature e.g. Brno nomenclature, Transcription Factor names etc.
@colour  #FFCC66

@column feature_type_id   Primary key, internal ID
@column analysis_id       @link analysis table ID
@column name              Name of feature_type
@column class             Class of feature_type
@column description       Text description
@column so_accession      Sequence ontology accession
@column so_name           Sequence ontology name
@column production_name   Name used in production

@see analysis
*/

DROP TABLE IF EXISTS `feature_type`;
CREATE TABLE `feature_type` (
  `feature_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(40) NOT NULL,
  `class` enum('Insulator','DNA','Regulatory Feature','Histone','RNA','Polymerase','Transcription Factor','Transcription Factor Complex','Regulatory Motif','Enhancer','Expression','Pseudo','Open Chromatin','Search Region','Association Locus','Segmentation State','DNA Modification','Transcription Start Site') DEFAULT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `so_accession` varchar(64) DEFAULT NULL,
  `so_name` varchar(255) DEFAULT NULL,
  `production_name` VARCHAR(120) DEFAULT NULL,
  PRIMARY KEY (`feature_type_id`),
  UNIQUE KEY `name_class_analysis_idx` (`name`,`class`,`analysis_id`),
  KEY `so_accession_idx` (`so_accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  associated_feature_type
@desc   Link table providing many to many mapping for feature_type entries.
@colour  #FFCC66

@column table_id         Internal table_id of linked table
@column feature_type_id  Internal table_id of linked @link feature_type
@column table_name       Name of linked table

@see  feature_type
*/

DROP TABLE IF EXISTS `associated_feature_type`;
CREATE TABLE `associated_feature_type` (
  `table_id` int(10) unsigned NOT NULL,
  `table_name` enum('annotated_feature','external_feature','regulatory_feature','feature_type') NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`table_id`,`table_name`,`feature_type_id`),
  KEY `feature_type_index` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**

@header  Set tables
@desc    Sets are containers for distinct sets of raw and/or processed data.
@colour  #66CCFF
@legend  #66CCFF Set tables
*/

/**
@table  feature_set
@desc   Container for genomic features defined by the result of an analysis e.g. peaks calls or regulatory features.
@colour  #66CCFF

@column feature_set_id  Internal ID
@column analysis_id     @link analysis ID
@column feature_type_id @link feature_type ID
@column name            Name for this feature set
@column type            Type of features contained e.g. external
@column description     Text description
@column display_label   Shorter more readable version of name

@see analysis
@see feature_type
*/

DROP TABLE IF EXISTS `feature_set`;
CREATE TABLE `feature_set` (
  `feature_set_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `feature_type_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `name` varchar(100) DEFAULT NULL,
  `type` enum('annotated','regulatory','external','segmentation','mirna_target') DEFAULT NULL,
  `description` varchar(80) DEFAULT NULL,
  `display_label` varchar(80) DEFAULT NULL,
  PRIMARY KEY (`feature_set_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `feature_type_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  peak_calling_qc_prop_reads_in_peaks
@desc
@colour  #66CCFF

@column peak_calling_qc_prop_reads_in_peaks_id  Internal ID
@column analysis_id                            @link analysis ID
@column peak_calling_id                         @link peak_calling ID
@column prop_reads_in_peaks                    Proportion of reads in peaks
@column total_reads                            Total number of reads
@column path
@column bam_file

@see analysis
@see peak_calling
*/

DROP TABLE IF EXISTS `peak_calling_qc_prop_reads_in_peaks`;
CREATE TABLE `peak_calling_qc_prop_reads_in_peaks` (
  `peak_calling_qc_prop_reads_in_peaks_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `peak_calling_id` int(10) unsigned NOT NULL,
  `prop_reads_in_peaks` double DEFAULT NULL,
  `total_reads` int(10) DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  `bam_file` varchar(512) NOT NULL,
  PRIMARY KEY (`peak_calling_qc_prop_reads_in_peaks_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  alignment

@desc   Alignment of reads from a ChIP-seq or similar experiment

@colour  #66CCFF

@column alignment_id     Internal ID
@column analysis_id      @link analysis ID The aligner used to create this alignment.
@column name             Name of the alignment.
@column bam_file_id      This is the data_file_id in the @link data_file for the bam file of this alignment.
@column bigwig_file_id   This is the data_file_id in the @link data_file for the bigwig file of this alignment.

@see analysis
@see data_file

*/

DROP TABLE IF EXISTS `alignment`;
CREATE TABLE `alignment` (
  `alignment_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `name` varchar(100) DEFAULT NULL,
  `bam_file_id` int(11) DEFAULT NULL,
  `bigwig_file_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`alignment_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  alignment_read_file

@desc   Linking table to connect alignments to the reads that were aligned.

@colour  #66CCFF

@column alignment_read_file_id     Internal ID
@column alignment_id     @link alignment ID
@column read_file_id     @link read_file ID

@see alignment
@see read_file

*/
DROP TABLE IF EXISTS `alignment_read_file`;
CREATE TABLE `alignment_read_file` (
  `alignment_read_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `alignment_id` int(10) unsigned NOT NULL,
  `read_file_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`alignment_read_file_id`,`alignment_id`),
  UNIQUE KEY `rset_table_idname_idx` (`alignment_id`,`read_file_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  alignment_qc_chance
@desc
@colour  #66CCFF

@column alignment_qc_chance        Internal ID
@column signal_result_set_id       @link alignment ID
@column analysis_id                @link analysis ID
@column p
@column q
@column divergence
@column z_score
@column percent_genome_enriched
@column input_scaling_factor
@column differential_percentage_enrichment
@column control_enrichment_stronger_than_chip_at_bin
@column first_nonzero_bin_at
@column pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome
@column path

@see alignment
@see analysis
*/

DROP TABLE IF EXISTS `alignment_qc_chance`;
CREATE TABLE `alignment_qc_chance` (
  `alignment_qc_chance` int(10) unsigned NOT NULL,
  `signal_alignment_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `p` double DEFAULT NULL,
  `q` double DEFAULT NULL,
  `divergence` double DEFAULT NULL,
  `z_score` double DEFAULT NULL,
  `percent_genome_enriched` double DEFAULT NULL,
  `input_scaling_factor` double DEFAULT NULL,
  `differential_percentage_enrichment` double DEFAULT NULL,
  `control_enrichment_stronger_than_chip_at_bin` double DEFAULT NULL,
  `first_nonzero_bin_at` double DEFAULT NULL,
  `pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome` double DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  PRIMARY KEY (`alignment_qc_chance`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  alignment_qc_flagstats
@desc
@colour  #66CCFF

@column alignment_qc_flagstats_id  Internal ID
@column alignment_id               @link alignment ID
@column analysis_id                @link analysis ID
@column category
@column qc_passed_reads
@column qc_failed_reads
@column path
@column bam_file

@see alignment
@see analysis
*/

DROP TABLE IF EXISTS `alignment_qc_flagstats`;
CREATE TABLE `alignment_qc_flagstats` (
  `alignment_qc_flagstats_id` int(10) unsigned NOT NULL,
  `alignment_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `category` varchar(100) NOT NULL,
  `qc_passed_reads` int(10) unsigned DEFAULT NULL,
  `qc_failed_reads` int(10) unsigned DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  `bam_file` varchar(512) NOT NULL,
  PRIMARY KEY (`alignment_qc_flagstats_id`),
  UNIQUE KEY `name_exp_idx` (`alignment_qc_flagstats_id`,`category`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  alignment_qc_phantom_peak
@desc
@colour  #66CCFF

@column alignment_qc_phantom_peak_id  Internal ID
@column analysis_id        @link analysis ID
@column alignment_id       @link alignment ID
@column filename
@column numReads
@column estFragLen
@column estFragLen2
@column estFragLen3
@column corr_estFragLen
@column corr_estFragLen2
@column corr_estFragLen3
@column phantomPeak
@column corr_phantomPeak
@column argmin_corr
@column min_corr
@column NSC
@column RSC
@column QualityTag
@column path

@see analysis
@see alignment
*/

DROP TABLE IF EXISTS `alignment_qc_phantom_peak`;
CREATE TABLE `alignment_qc_phantom_peak` (
  `alignment_qc_phantom_peak_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `alignment_id` int(10) unsigned NOT NULL,
  `filename` varchar(512) NOT NULL,
  `numReads` int(10) unsigned NOT NULL,
  `estFragLen` double DEFAULT NULL,
  `estFragLen2` double DEFAULT NULL,
  `estFragLen3` double DEFAULT NULL,
  `corr_estFragLen` double DEFAULT NULL,
  `corr_estFragLen2` double DEFAULT NULL,
  `corr_estFragLen3` double DEFAULT NULL,
  `phantomPeak` int(10) unsigned NOT NULL,
  `corr_phantomPeak` double DEFAULT NULL,
  `argmin_corr` int(10) DEFAULT NULL,
  `min_corr` double DEFAULT NULL,
  `NSC` double DEFAULT NULL,
  `RSC` double DEFAULT NULL,
  `QualityTag` int(10) DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  PRIMARY KEY (`alignment_qc_phantom_peak_id`),
  KEY `filename_idx` (`filename`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  data_file
@colour  #66CCFF
@desc

@column data_file_id   Internal ID
@column table_id       Primary key of linked entity e.g. @link external_feature_file or @link segmentation_file or @link alignment
@column table_name     Name of linked table (external_feature_file, segmentation_file, alignment)
@column path           Either a full filepath or a directory which the API will use to build the filepath
@column file_type      Type of data file ('BAM','BAMCOV','BIGBED','BIGWIG','VCF','CRAM','DIR')
@column md5sum         md5sum of data file

@see external_feature_file
@see segmentation_file
@see alignment
*/


DROP TABLE IF EXISTS `data_file`;
CREATE TABLE `data_file` (
  `data_file_id` int(11) NOT NULL AUTO_INCREMENT,
  `table_id` int(10) unsigned NOT NULL,
  `table_name` varchar(32) NOT NULL,
  `path` varchar(255) NOT NULL,
  `file_type` enum('BAM','BAMCOV','BIGBED','BIGWIG','VCF','CRAM','DIR') NOT NULL DEFAULT 'BAM',
  `md5sum` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`data_file_id`),
  UNIQUE KEY `table_id_name_path_idx` (`table_id`,`table_name`,`path`),
  UNIQUE KEY `data_file_id` (`data_file_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  read_file
@desc
@colour  #66CCFF

@column read_file_id    Internal ID
@column name            Name for the read file object
@column analysis_id     @link analysis ID
@column is_paired_end   Indicates whether it is paired end
@column file_size       (not used)
@column read_length     (not used)
@column md5sum          (not used)
@column file            Location of the read file on disk (not public)
@column notes           (not used)

@see read_file_experimental_configuration
@see analysis

*/

DROP TABLE IF EXISTS `read_file`;
CREATE TABLE `read_file` (
  `read_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(300) NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `is_paired_end` tinyint(1) DEFAULT NULL,
  `file_size` bigint(20) DEFAULT NULL,
  `read_length` int(10) DEFAULT NULL,
  `md5sum` varchar(45) DEFAULT NULL,
  `file` text,
  `notes` text,
  PRIMARY KEY (`read_file_id`),
  UNIQUE KEY `read_file_id_idx` (`read_file_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  read_file_experimental_configuration
@desc
@colour  #66CCFF

@column read_file_experimental_configuration_id  Internal ID
@column read_file_id                             @link read_file id of the read file that is being described.
@column experiment_id                            @link experiment id of the experiment during which the read file was generated.
@column biological_replicate                     @link the biological replicate number
@column technical_replicate                      @link the technical replicate number
@column paired_end_tag                           (not used yet)
@column multiple                                 (not used yet)

@see read_file
@see experiment

*/
DROP TABLE IF EXISTS `read_file_experimental_configuration`;
CREATE TABLE `read_file_experimental_configuration` (
  `read_file_experimental_configuration_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `read_file_id` int(10) unsigned DEFAULT NULL,
  `experiment_id` int(10) unsigned NOT NULL,
  `biological_replicate` tinyint(3) unsigned NOT NULL DEFAULT '1',
  `technical_replicate` tinyint(3) unsigned NOT NULL DEFAULT '1',
  `paired_end_tag` int(11) DEFAULT NULL,
  `multiple` int(11) DEFAULT '1',
  PRIMARY KEY (`read_file_experimental_configuration_id`),
  UNIQUE KEY `name_exp_idx` (`experiment_id`,`biological_replicate`,`technical_replicate`),
  KEY `experiment_idx` (`experiment_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@header Array design tables
@colour  #FF6666
@legend  #FF6666   Array design tables
*/


/**
@table  array
@desc   Contains information defining an array or array set.
@colour  #FF6666

@column array_id    Internal ID
@column name        Name of array
@column format      Format of array e.g. EXPRESSION, TILED
@column vendor      Name of array vendor e.g. AFFY
@column description Text description
@column type        Array type e.g. OLIGO, PCR
@column class       Array class e.g. AFFY_ST, ILLUMINA_INFINIUM

@column is_probeset_array         Indicates whether the array is organised into probe sets.
@column is_linked_array
@column has_sense_interrogation   Indicates whether the array has sense interrogation

*/

DROP TABLE IF EXISTS `array`;
CREATE TABLE `array` (
  `array_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(40) DEFAULT NULL,
  `format` varchar(20) DEFAULT NULL,
  `vendor` varchar(40) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `type` varchar(20) DEFAULT NULL,
  `class` varchar(20) DEFAULT NULL,
  `is_probeset_array` tinyint(1) NOT NULL DEFAULT '0',
  `is_linked_array` tinyint(1) NOT NULL DEFAULT '0',
  `has_sense_interrogation` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`array_id`),
  UNIQUE KEY `vendor_name_idx` (`vendor`,`name`),
  UNIQUE KEY `class_name_idx` (`class`,`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  array_chip
@desc   Represents the individual array chip design as part of an array or array set.
@colour  #FF6666

@column array_chip_id   Internal ID
@column array_id        @link array ID
@column design_id       ID/Accession defined by vendor
@column name            Name of array_chip

@see array
*/

DROP TABLE IF EXISTS `array_chip`;
CREATE TABLE `array_chip` (
  `array_chip_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `design_id` varchar(100) DEFAULT NULL,
  `array_id` int(10) unsigned NOT NULL,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`array_chip_id`),
  UNIQUE KEY `array_design_idx` (`array_id`,`design_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  probe_set
@desc   The table contains information about probe sets.
@colour  #FF6666

@column probe_set_id Internal ID
@column name          Name of the probe set
@column size          Integer size of the probe set i.e. how many probe is contains
@column family        Generic descriptor for probe_set e.g. ENCODE_REGIONS, RANDOM etc. Not used
@column array_chip_id @link array_chip ID of the array chip to which this probe set belongs.

@see array_chip

*/

DROP TABLE IF EXISTS `probe_set`;
CREATE TABLE `probe_set` (
  `probe_set_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) NOT NULL,
  `size` smallint(6) unsigned NOT NULL,
  `family` varchar(20) DEFAULT NULL,
  `array_chip_id` int(10) DEFAULT NULL,
  PRIMARY KEY (`probe_set_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  probe_set_transcript
@desc   This table maps probe sets to transcripts.
@colour  #FF6666

@column probe_set_transcript_id Internal ID
@column probe_set_id            Id of the @link probe_set
@column stable_id               Stable id of the transcript to which it has been mapped
@column description             Details about the mapping as text

@see probe_set

*/

DROP TABLE IF EXISTS `probe_set_transcript`;
CREATE TABLE `probe_set_transcript` (
  `probe_set_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_set_id` int(10) unsigned NOT NULL,
  `stable_id` varchar(128) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_set_transcript_id`),
  KEY `probe_set_transcript_id_idx` (`probe_set_transcript_id`),
  KEY `probe_set_transcript_stable_id_idx` (`stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  probe
@desc   Defines individual probe designs across one or more array_chips.
@colour  #FF6666

@column probe_id      Internal ID
@column array_chip_id @link array_chip ID
@column probe_set_id  @link probe_set ID
@column name          Name of the probe set
@column length        Integer bp length of the probe
@column class         Class of the probe e.g. CONTROL, EXPERIMENTAL etc.
@column description   Text description
@column probe_seq_id  @link probe_seq ID

@see array_chip
@see probe_set
@see probe_seq
*/

DROP TABLE IF EXISTS `probe`;
CREATE TABLE `probe` (
  `probe_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_set_id` int(10) unsigned DEFAULT NULL,
  `name` varchar(100) NOT NULL,
  `length` smallint(6) unsigned NOT NULL,
  `array_chip_id` int(10) unsigned NOT NULL,
  `class` varchar(20) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `probe_seq_id` int(10) DEFAULT NULL,
  PRIMARY KEY (`probe_id`,`name`,`array_chip_id`),
  KEY `probe_set_idx` (`probe_set_id`),
  KEY `array_chip_idx` (`array_chip_id`),
  KEY `name_idx` (`name`),
  KEY `probe_seq_idx` (`probe_seq_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  probe_seq
@desc   Probe sequences
@colour  #FF6666

@column probe_seq_id          Internal ID
@column sequence              Probe sequence
@column sequence_upper        Probe sequence in uppcase letters.
@column sequence_upper_sha1   Sha1 hashsum of the uppcase of the probe sequence.

@see probe
*/

DROP TABLE IF EXISTS `probe_seq`;
CREATE TABLE `probe_seq` (
  `probe_seq_id` int(10) NOT NULL AUTO_INCREMENT,
  `sequence` text NOT NULL,
  `sequence_upper` text NOT NULL,
  `sequence_upper_sha1` char(40) NOT NULL,
  PRIMARY KEY (`probe_seq_id`),
  UNIQUE KEY `sequence_upper_sha1` (`sequence_upper_sha1`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  probe_transcript
@desc   This table maps probes to transcripts.
@colour  #FF6666

@column probe_transcript_id Internal ID
@column probe_id            Id of the @link probe_set
@column stable_id           Stable id of the transcript to which it has been mapped
@column description         Details about the mapping as text

@see probe

*/

DROP TABLE IF EXISTS `probe_transcript`;
CREATE TABLE `probe_transcript` (
  `probe_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_id` int(10) unsigned NOT NULL,
  `stable_id` varchar(128) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_transcript_id`),
  KEY `probe_transcript_id` (`probe_transcript_id`),
  KEY `probe_transcript_stable_id_idx` (`stable_id`),
  KEY `probe_transcript_idx` (`probe_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@header Experiment tables
@desc   These define the experimental meta and raw data.
@colour  #00FF00
@legend  #00FF00 Experiment tables
*/


/**
@table  experiment
@desc   Represents a sequencing experiment. Sequencing runs (input_subsets) link to this.
@colour  #00FF00

@column experiment_id           Internal ID
@column name                    Name of experiment
@column experimental_group_id   @link experimental_group ID
@column control_id              @link experiment ID
@column is_control              Boolean, true means that this experiment is a control.
@column feature_type_id         @link feature_type table ID
@column epigenome_id            @link epigenome ID
@column archive_id              ENA experiment identifier enabling access to specific raw data

@see epigenome
@see experimental_group
@see feature_type
*/

DROP TABLE IF EXISTS `experiment`;
CREATE TABLE `experiment` (
  `experiment_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  `experimental_group_id` smallint(6) unsigned DEFAULT NULL,
  `control_id` int(10) unsigned DEFAULT NULL,
  `is_control` tinyint(3) unsigned DEFAULT '0',
  `feature_type_id` int(10) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `archive_id` varchar(60) DEFAULT NULL,
  PRIMARY KEY (`experiment_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `experimental_group_idx` (`experimental_group_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `epigenome_idx` (`epigenome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  experimental_group
@desc   Think: Consortium or laboratory that produced sequencing experiments (@see experiment).
@colour  #00FF00

@column experimental_group_id  Internal ID
@column contact                Contact details e.g. email
@column description            Text description
@column is_project             Large or small scale project
@column location               Geographic location of group
@column name                   Name of group
@column url                    Url for Project page

*/

DROP TABLE IF EXISTS `experimental_group`;
CREATE TABLE `experimental_group` (
  `experimental_group_id` smallint(6) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(40) NOT NULL,
  `location` varchar(120) DEFAULT NULL,
  `contact` varchar(40) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `url` varchar(255) DEFAULT NULL,
  `is_project` tinyint(1) DEFAULT '0',
  PRIMARY KEY (`experimental_group_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@header  Ancilliary tables
@desc    These contain data types which are used across many of the above tables and
         are quite often denormalised to store generic associations to several table,
         this avoids the need for multiple sets of similar tables. Some of these tables
         have been omitted from the schema diagram.
@colour  #808000
*/


/**
@table  epigenome
@desc   The epigenomes known in Ensembl regulation.
@colour  #808000

@column epigenome_id         Internal ID
@column name                 Name of the epigenome
@column display_label        Name of epigenome for displaying on the website
@column description          Text description, used in the z-menu that appears when hovering over the epigenome name
@column production_name      Production name of the epigenome
@column gender               Gender i.e. 'male', 'female', 'hermaphrodite', 'unknown' or 'mixed'
@column ontology_accession   External accession id
@column ontology             The resource the ontology_accession refers to, currently either EFO or CL
@column tissue               Tissue origin/type

*/

DROP TABLE IF EXISTS `epigenome`;
CREATE TABLE `epigenome` (
  `epigenome_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(120) NOT NULL,
  `display_label` varchar(30) NOT NULL,
  `description` varchar(80) DEFAULT NULL,
  `production_name` varchar(120) DEFAULT NULL,
  `gender` enum('male','female','hermaphrodite','mixed','unknown') DEFAULT 'unknown',
  `ontology_accession` varchar(20) DEFAULT NULL,
  `ontology` enum('EFO','CL') DEFAULT NULL,
  `tissue` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`epigenome_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**

@header  Core tables
@desc    These are exact clones of the corresponding core schema tables, hence have been omitted
         from the schema diagram. See <a href='../core/core_schema.html'>core schema docs</a> for more details.
@colour  #000000
@legend  #000000 Core/Core like tables (Omited from schema diagram)
*/

/**
@table analysis
@desc Usually describes a program and some database that together are used to create a feature on a piece of sequence.
Each feature is marked with an analysis_id. The most important column is logic_name, which is used by the webteam to render a feature
 correctly on contigview (or even retrieve the right feature).
Logic_name is also used in the pipeline to identify the analysis which has to run in a given status of the pipeline.
The module column tells the pipeline which Perl module does the whole analysis, typically a RunnableDB module.
@colour  #000000

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
  `analysis_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `created` datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  `logic_name` varchar(100) NOT NULL,
  `db` varchar(120) DEFAULT NULL,
  `db_version` varchar(40) DEFAULT NULL,
  `db_file` varchar(120) DEFAULT NULL,
  `program` varchar(80) DEFAULT NULL,
  `program_version` varchar(40) DEFAULT NULL,
  `program_file` varchar(80) DEFAULT NULL,
  `parameters` text,
  `module` varchar(80) DEFAULT NULL,
  `module_version` varchar(40) DEFAULT NULL,
  `gff_source` varchar(40) DEFAULT NULL,
  `gff_feature` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`analysis_id`),
  UNIQUE KEY `logic_name_idx` (`logic_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table analysis_description
@desc Allows the storage of a textual description of the analysis, as well as a "display label", primarily for the EnsEMBL web site.
@colour  #000000

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
  `displayable` tinyint(1) NOT NULL DEFAULT '1',
  `web_data` text,
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table meta
@desc Stores data about the data in the current schema. Unlike other tables, data in the meta table is stored as key-value pairs. These data include details about the database, RegulatoryBuild and patches.  The species_id field of the meta table is used in multi-species databases and makes it possible to have species-specific meta key-value pairs.  The species-specific meta key-value pairs needs to be repeated for each species_id.  Entries in the meta table that are not specific to any one species, such as the schema.version key and any other schema-related information must have their species_id field set to NULL
. The default species_id, and the only species_id value allowed in single-species databases, is 1.
@colour  #000000


@column meta_id                    Internal identifier.
@column species_id                 Indentifies the species for multi-species databases.
@column meta_key                   Name of the meta entry, e.g. "schema_version".
@column meta_value                 Corresponding value of the key, e.g. "61".

*/


DROP TABLE IF EXISTS `meta`;
CREATE TABLE `meta` (
  `meta_id` int(10) NOT NULL AUTO_INCREMENT,
  `species_id` int(10) unsigned DEFAULT '1',
  `meta_key` varchar(46) NOT NULL,
  `meta_value` varchar(950) NOT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `species_key_value_idx` (`species_id`,`meta_key`,`meta_value`),
  KEY `species_value_idx` (`species_id`,`meta_value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- Add necessary meta values
INSERT INTO meta (meta_key, meta_value) VALUES ('schema_type', 'funcgen');

-- Update and remove these for each release to avoid erroneous patching
INSERT INTO `meta` (species_id, meta_key, meta_value) VALUES (null, 'schema_version', '92');
INSERT INTO `meta` (species_id, meta_key, meta_value) VALUES (null, 'patch','patch_91_92_a.sql|schema_version');
INSERT INTO `meta` (species_id, meta_key, meta_value) VALUES (null, 'patch','patch_91_92_b.sql|Drop column paired_with from table read_file');
INSERT INTO `meta` (species_id, meta_key, meta_value) VALUES (NULL, 'patch','patch_91_92_c.sql|Create underlying_structure table');

/**
@table meta_coord
@desc Describes which co-ordinate systems the different feature tables use.
@colour  #000000

@column table_name              Ensembl database table name.
@column coord_system_id         Table ID for @link coord_system
@column max_length              Longest sequence length.

*/

DROP TABLE IF EXISTS `meta_coord`;
CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `max_length` int(11) DEFAULT NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table  associated_xref
@desc   This table associates extra associated annotations with a given ontology xref evidence and source under a specific condition.   For GO this allows qualifiers (with/from) or annotation extensions to be added to a given ontology annotation.
@colour  #000000

@column associated_xref_id Associated xref id. Primary key, internal identifier
@column object_xref_id Object xref id this associated xref is linked to. Foreign key linked to the @link object_xref table
@column xref_id Xref which is the associated term. Foreign key linked to the @link xref table
@column source_xref_id Xref which is source of this association. Foreign key linked to the @link xref table
@column condition_type The type of condition this link occurs in e.g. evidence, from, residue or assigned_by
@column associated_group_id Foreign key to allow for @link associated_group
@column rank The rank in which the association occurs within an @link associated_group

@see object_xref
@see associcated_group
@see xref
*/
DROP TABLE IF EXISTS `associated_xref`;
CREATE TABLE `associated_xref` (
  `associated_xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `object_xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `source_xref_id` int(10) unsigned DEFAULT NULL,
  `condition_type` varchar(128) DEFAULT NULL,
  `associated_group_id` int(10) unsigned DEFAULT NULL,
  `rank` int(10) unsigned DEFAULT '0',
  PRIMARY KEY (`associated_xref_id`),
  UNIQUE KEY `object_associated_source_type_idx` (`object_xref_id`,`xref_id`,`source_xref_id`,`condition_type`,`associated_group_id`),
  KEY `associated_source_idx` (`source_xref_id`),
  KEY `associated_object_idx` (`object_xref_id`),
  KEY `associated_idx` (`xref_id`),
  KEY `associated_group_idx` (`associated_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table  associated_group
@desc   Groups together xref associations under a single description. Used when more than one associated xref term must be used to describe a condition
@colour  #000000

@column associated_group_id Associated group id. Primary key, internal identifier
@column description Optional description for this group
*/
DROP TABLE IF EXISTS `associated_group`;
CREATE TABLE `associated_group` (
  `associated_group_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `description` varchar(128) DEFAULT NULL,
  PRIMARY KEY (`associated_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table identity_xref
@desc Describes how well a particular xref object matches the EnsEMBL object.
@colour  #000000

@column object_xref_id     Foreign key references to the @link object_xref table.
@column xref_identity      Percentage identity.
@column ensembl_identity   Percentage identity.
@column xref_start         Xref sequence start.
@column xref_end           Xref sequence end.
@column ensembl_start      Ensembl sequence start.
@column ensembl_end        Ensembl sequence end.
@column cigar_line         Used to encode gapped alignments.
@column score              Match score.
@column evalue             Match evalue.

@see object_xref

*/


DROP TABLE IF EXISTS `identity_xref`;
CREATE TABLE `identity_xref` (
  `object_xref_id` int(10) unsigned NOT NULL,
  `xref_identity` int(5) DEFAULT NULL,
  `ensembl_identity` int(5) DEFAULT NULL,
  `xref_start` int(11) DEFAULT NULL,
  `xref_end` int(11) DEFAULT NULL,
  `ensembl_start` int(11) DEFAULT NULL,
  `ensembl_end` int(11) DEFAULT NULL,
  `cigar_line` text,
  `score` double DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  PRIMARY KEY (`object_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table external_synonym
@desc Some xref objects can be referred to by more than one name. This table relates names to xref IDs.
@colour  #000000

@column xref_id           Foreign key references @link xref table
@column synonym           Synonym

@see xref

*/


DROP TABLE IF EXISTS `external_synonym`;
CREATE TABLE `external_synonym` (
  `xref_id` int(10) unsigned NOT NULL,
  `synonym` varchar(100) NOT NULL,
  PRIMARY KEY (`xref_id`,`synonym`),
  KEY `name_index` (`synonym`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=20;


/**
@table external_db
@desc Stores data about the external databases in which the objects described in the xref table are stored.
@colour  #000000

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


DROP TABLE IF EXISTS `external_db`;
CREATE TABLE `external_db` (
  `external_db_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `db_name` varchar(100) NOT NULL,
  `db_release` varchar(255) DEFAULT NULL,
  `status` enum('KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO') NOT NULL,
  `dbprimary_acc_linkable` tinyint(1) NOT NULL DEFAULT '1',
  `priority` int(11) NOT NULL,
  `db_display_name` varchar(255) DEFAULT NULL,
  `type` enum('ARRAY','ALT_TRANS','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL') DEFAULT NULL,
  `secondary_db_name` varchar(255) DEFAULT NULL,
  `secondary_db_table` varchar(255) DEFAULT NULL,
  `description` text,
  PRIMARY KEY (`external_db_id`),
  UNIQUE KEY `db_name_release_idx` (`db_name`,`db_release`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=80;


/**
@table ontology_xref
@desc This table associates ontology terms/accessions to Ensembl objects (primarily EFO/SO). NOTE: Currently not in use
@colour  #000000

@column object_xref_id          Foreign key references to the @link object_xref table.
@column source_xref_id          Foreign key references to the @link xref table.
@column linkage_type            Defines type of linkage

@see object_xref

*/


DROP TABLE IF EXISTS `ontology_xref`;
CREATE TABLE `ontology_xref` (
  `object_xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `source_xref_id` int(10) unsigned DEFAULT NULL,
  `linkage_type` enum('IC','IDA','IEA','IEP','IGI','IMP','IPI','ISS','NAS','ND','TAS','NR','RCA') NOT NULL,
  UNIQUE KEY `object_xref_id_2` (`object_xref_id`,`source_xref_id`,`linkage_type`),
  KEY `object_xref_id` (`object_xref_id`),
  KEY `source_xref_id` (`source_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- This is just an empty to table to avoid having to rework all the core sql and API to accomodate eFG specific xref schema

/**
@table unmapped_reason
@desc Describes the reason why a mapping failed.
@colour  #000000

@column unmapped_reason_id           Internal identifier.
@column summary_description          Summarised description.
@column full_description             Full description.

@see unmapped_object
*/

DROP TABLE IF EXISTS `unmapped_reason`;
CREATE TABLE `unmapped_reason` (
  `unmapped_reason_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `summary_description` varchar(255) DEFAULT NULL,
  `full_description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`unmapped_reason_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**

@header  Core like tables
@desc    These are almost exact clones of the corresponding core schema tables. Some contain extra fields or different enum values to support the funcgen schema. These have been omitted from the schema diagram.
@colour  #000000

*/

/**
@table xref
@desc Holds data about objects which are external to EnsEMBL, but need to be associated with EnsEMBL objects.
Information about the database that the external object is stored in is held in the external_db table entry referred to by the external_db column.

@column xref_id                 Primary key, internal identifier.
@column external_db_id          Foreign key references to the @link external_db table.
@column dbprimary_acc           Primary accession number.
@column display_label           Display label for the EnsEMBL web site.
@column version                 Object version.
@column description             Object description.
@column info_type               'PROJECTION', 'MISC', 'DEPENDENT','DIRECT', 'SEQUENCE_MATCH','INFERRED_PAIR', 'PROBE','UNMAPPED', 'COORDINATE_OVERLAP', 'CHECKSUM'.
@column info_text               Text

@see external_db
@see external_synonym

*/


DROP TABLE IF EXISTS `xref`;
CREATE TABLE `xref` (
  `xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `external_db_id` int(10) unsigned DEFAULT NULL,
  `dbprimary_acc` varchar(512) NOT NULL,
  `display_label` varchar(512) NOT NULL,
  `version` varchar(10) DEFAULT NULL,
  `description` text,
  `info_type` enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','COORDINATE_OVERLAP','CHECKSUM') NOT NULL DEFAULT 'NONE',
  `info_text` varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (`xref_id`),
  UNIQUE KEY `id_index` (`dbprimary_acc`,`external_db_id`,`info_type`,`info_text`,`version`),
  KEY `display_index` (`display_label`),
  KEY `info_type_idx` (`info_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=100;



/**
@table object_xref
@desc Describes links between Ensembl objects and objects held in external databases.  The Ensembl object can be one of several types; the type is held in the ensembl_object_type column.  The ID of the particular Ensembl gene, translation or whatever is given in the ensembl_id column.  The xref_id points to the entry in the xref table that holds data about the external object.  Each Ensembl object can be associated with zero or more xrefs. An xref object can be associated with one or more Ensembl objects.
@colour  #000000

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
CREATE TABLE `object_xref` (
  `object_xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ensembl_id` int(10) unsigned NOT NULL,
  `ensembl_object_type` enum('Epigenome','Experiment','RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType','MirnaTargetFeature','ProbeSet','Probe','ProbeFeature') NOT NULL,
  `xref_id` int(10) unsigned NOT NULL,
  `linkage_annotation` varchar(255) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  PRIMARY KEY (`object_xref_id`),
  UNIQUE KEY `xref_idx` (`xref_id`,`ensembl_object_type`,`ensembl_id`,`analysis_id`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `ensembl_idx` (`ensembl_object_type`,`ensembl_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=40;

/**
@table unmapped_object
@desc Describes why a particular external entity was not mapped to an ensembl one.
@colour  #000000

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
  `unmapped_object_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('xref','probe2transcript','array_mapping') NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `external_db_id` int(10) unsigned DEFAULT NULL,
  `identifier` varchar(255) NOT NULL,
  `unmapped_reason_id` int(10) unsigned NOT NULL,
  `query_score` double DEFAULT NULL,
  `target_score` double DEFAULT NULL,
  `ensembl_id` int(10) unsigned DEFAULT '0',
  `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType','Probe','ProbeSet','ProbeFeature') NOT NULL,
  `parent` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`unmapped_object_id`),
  UNIQUE KEY `unique_unmapped_obj_idx` (`ensembl_id`,`ensembl_object_type`,`identifier`,`unmapped_reason_id`,`parent`,`external_db_id`),
  KEY `anal_exdb_idx` (`analysis_id`,`external_db_id`),
  KEY `id_idx` (`identifier`(50)),
  KEY `ext_db_identifier_idx` (`external_db_id`,`identifier`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/**
@table underlying_structure
@desc Associates regulatory features to motif features
@colour  #000000

@column underlying_structure_id    Internal identifier.
@column regulatory_feature_id      Foreign key references to the @link regulatory_feature table.
@column motif_feature_id           Foreign key references to the @link motif_feature table.

@see regulatory_feature
@see underlying_structure
*/

DROP TABLE IF EXISTS `underlying_structure`;
CREATE TABLE `underlying_structure` (
  `underlying_structure_id` int(11) NOT NULL AUTO_INCREMENT,
  `regulatory_feature_id` int(11) NOT NULL,
  `motif_feature_id` int(11) NOT NULL,
  PRIMARY KEY (`underlying_structure_id`)
)ENGINE=MyISAM DEFAULT CHARSET=latin1;