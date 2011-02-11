/** 
@header patch_61_62_e.sql - dbfile_registry
@desc   Add generic dbfile_registry table to store filepaths

@table  dbfile_registry
@desc   This generic table contains a simple registry of paths to support 
        flat file(DBFile) access. This should be left joined from the
        relevant adaptor e.g. ResultSetAdaptor
@column table_id   Primary key of linked dbfile entity e.g. result_set or 
        analysis
@column table_name Name of linked table
@column path       Either a full filepath or a directory which the API will 
                   use to build the filepath
*/


DROP TABLE IF EXISTS `dbfile_registry`;
CREATE TABLE `dbfile_registry` (
   `table_id` int(10) unsigned NOT NULL,
   `table_name` varchar(32)NOT NULL,	
   `path` varchar(255) NOT NULL,
   PRIMARY KEY  (`table_id`, `table_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;



# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_61_62_e.sql|dbfile_registry');


