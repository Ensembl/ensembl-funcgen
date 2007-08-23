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
delete from supporting_set where member_set_id=0;


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
  `attribute_feature_type` enum('annotated', 'supporting') default NULL,
  PRIMARY KEY  (`regulatory_feature_id`, `attribute_feature_type`, `attribute_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=17;


-- do we need key on feature_type and or feature_id?

-- change all ids to int(10), 
-- add UNIQUE KEY `name_idx` (name) on data_set
-- alter key on array_chip  UNIQUE KEY `array_design_idx` (`array_id`, `design_id`)
-- recreate predicted_feature as annotated_feature

-- add enum on channel type TOTAL, EXPERIMENTAL & psuedo channels

-- add core regulatory tables as supporting_feature tables
-- regulatory_factor_coding is empty and unused?


-- supporting feature table or import directly into annotated feature
-- where do we draw the line between what goes in supporting rather than annotated?
-- what do we do about the overloading of the feature_type table?
-- supporting features must have multiple insatnces of feature_type with unique ids
-- i.e. high volume e.g. individual miRNAs

CREATE TABLE `supporting_feature` (
 `supporting_feature_id` int(10) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `seq_region_id` int(10) unsigned NOT NULL default '0',
  `seq_region_start` int(10) unsigned NOT NULL default '0',
  `seq_region_end` int(10) unsigned NOT NULL default '0',
  `seq_region_strand` tinyint(4) NOT NULL default '0',
  `analysis_id` smallint(5) unsigned NOT NULL default '0',
  `regulatory_factor_id` int(10) unsigned default NULL,
  `coord_system_id` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`supporting_feature_id`),
  KEY `seq_region_idx` (`seq_region_id`,`analysis_id`,`seq_region_start`),
  KEY `seq_region_idx_2` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- imported core regulatory_feature/factor/search_region here
--insert into supporting_feature(select *, 1 from regulatory_feature);

--strip off ENST prefixes for names??


-- add feature_type_id column to replace type
alter table regulatory_factor add feature_type_id int(10) unsigned NOT NULL default '0';
insert into feature_type(name, class, description) values('miRNA Target', 'RNA', 'miRNA target motif');
insert into feature_type(name, class, description) values('Transcription Factor', 'TRANSCRIPTION FACTOR', 'Transcription factor motif');
insert into feature_type(name, class, description) values('Transcription Factor Complex', 'TRANSCRIPTION FACTOR', 'Transcription complex factor motif');

update regulatory_factor rf, feature_type ft set rf.feature_type_id=ft.feature_type_id where rf.type='miRNA_target' and ft.name='miRNA Target';

--don't need to update other as they are all NULL

-- remove old type column

alter table regulatory_factor drop type;


--we need to split the regulatory_search_region table to extract the xrefs


-- add regulatory factor types to feature_type


-- need to finish off the reg feature stuff, but doing cs stuff first



