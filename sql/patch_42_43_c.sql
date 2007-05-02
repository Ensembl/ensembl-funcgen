-- experimental_chip drop description, add replicate, cell_type, feature_type
-- will need to hard code unique_ids to update cell_type_id and feature_type_id
-- cell_types add manuaully as they are species dependent

CREATE TABLE `cell_type` (
   `cell_type_id` int(10) unsigned NOT NULL auto_increment,
   `name`  varchar(120) not NULL,
   `display_label` varchar(20) default NULL,
   `description` varchar(40) default NULL,
   PRIMARY KEY  (`cell_type_id`),
   UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


-- add cell type dependent on species

alter table experimental_chip drop description;
alter table experimental_chip add feature_type_id int(10) unsigned default NULL;
alter table experimental_chip add cell_type_id int(10) unsigned default NULL;
alter table experimental_chip add replicate varchar(20) default 'UNKNOWN';


--update feature type manually? could do with experiment feature?
--update experimental_chip ec, experiment e set ec.cell_type_id=1 where e.name like "%HeLa" and e.experiment_id=ec.experiment_id;
--update experimental_chip ec, experiment e set ec.cell_type_id=2 where e.name like "%GM0%" and e.experiment_id=ec.experiment_id;

--update experimental_chip ec, experiment e set ec.feature_type_id=2 where e.name like "%H4ac%" and ec.experiment_id = e.experiment_id;
