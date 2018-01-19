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
@header patch_65_66_b.sql - cell_type.tissue and lineage
@desc   Add lineage table and tissue field to cell_type
*/

ALTER table cell_type ADD `tissue` varchar(50) default NULL;
analyze table cell_type;
optimize table cell_type;



/** 
@table  cell_type_lineage
@desc   Links cell_types to lineage terms
@colour  #808000

@column cell_type_id	Internal ID
@column lineage_id      Internal ID
@column most_specific   Denotes most specific term for this cell_type

@see cell_type
@see lineage
*/


DROP TABLE IF EXISTS `cell_type_lineage`;
CREATE TABLE `cell_type_lineage` (
   `cell_type_id` int(10) unsigned NOT NULL,
   `lineage_id` int(10) unsigned NOT NULL,
   `most_specific` boolean default NULL,
   PRIMARY KEY  (`cell_type_id`, `lineage_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- most_specific could be infered from lineage chain
-- add description?
-- most_specific here as this may be dependant on the cell_type


/** 
@table  lineage
@desc   Contains cell lineage information
@colour  #808000

@column cell_type_id	  Internal ID
@column name              Lineage name
@column efo_id            Experimental Factor Ontology ID
@column parent_lineage_id Internal ID of immediate parent term

@see cell_type_lineage
@see cell_type
*/


DROP TABLE IF EXISTS `lineage`;
CREATE TABLE `lineage` (
   `lineage_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(100) not NULL,
   `efo_id` varchar(20) DEFAULT NULL,
   `parent_lineage_id` int(10) unsigned DEFAULT NULL,
   PRIMARY KEY  (`lineage_id`),
   UNIQUE KEY `name_idx` (`name`),
   UNIQUE KEY `efo_idx` (`efo_id`),
   KEY `parent_linage_idx`(`parent_lineage_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- EFO and CL ontology (and hence ensembl_ontology DB) do not handle lineage
-- very well, so model here for now.
-- If parent_lineage_id is 0, then is effectively root lineage term.


-- Need to add xrefs to ensembl_ontology DB. Write HC for this eventually.


-- Add to foreign keys file!


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_65_66_b.sql|cell_type.tissue_and_lineage');


