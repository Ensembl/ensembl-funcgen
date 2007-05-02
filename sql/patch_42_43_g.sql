-- Add status_name table with all states present, alter status table accordingly

CREATE TABLE `status_name` (
   `status_name_id` int(10) unsigned NOT NULL auto_increment,
   `name` varchar(20) default NULL,	
   PRIMARY KEY  (`status_name_id`),
   UNIQUE KEY `status_name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


INSERT into status_name(name) SELECT distinct(state) from status;

ALTER table status ADD COLUMN status_name_id int(10);

UPDATE status s, status_name sn SET s.status_name_id=sn.status_name_id WHERE s.state=sn.name;

ALTER table status CHANGE status_name_id status_name_id int(10) NOT NULL;
ALTER table status DROP PRIMARY KEY;
ALTER table status ADD PRIMARY KEY (table_id, table_name, status_name_id);
ALTER table status DROP state;


-- now remove odl tables and add experimental_design_type link table
drop table experiment_prediction;
drop table experiment_feature_type;



