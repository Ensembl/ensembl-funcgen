# Foreign key relationships in the Ensembl schema (see efg.sql)
#
# This file is intended as a reference since some of the relationships
# are not obvious.
#
# Note that these constraints are not actually used by Ensembl for 
# performance reasons, and referential integrity is enforced at the
# application level. Also MySQL currently does not support foreign
# key constraints on MyISAM tables.


-- To run this script you must first create an empty DB and update all tables to InnoDB
-- mysqlw -pXXXXXXX -e 'DROP DATABASE innodb_funcgen_65'
-- mysqlw -pXXXXXXX -e 'CREATE DATABASE innodb_funcgen_65'
-- mysqlw -pXXXXXXX innodb_funcgen_65 < efg.sql
-- mysqlro innodb_funcgen_65 -N -e "show tables" | while read t; do if [[ $t != meta ]]; then echo "Altering $t to InnoDB"; mysqlw -pXXXXXXX innodb_funcgen_65 -N -e "ALTER TABLE $t engine=InnoDB"; fi; done

-- Then, optionally disable some to remove complexity when generating ER diagrams, before running:
-- mysqlw -pXXXXXXX innodb_funcgen_65 < foreign_keys.sql

-- Some potential errors
-- ERROR 1071 (42000) at line 1: Specified key was too long; max key length is 767 bytes
-- This is caused by meta table, which has no foreign keys anyway, so skip this.




-- Last updated for v65



### eFG FKs

-- segmentation_feature
ALTER table segmentation_feature ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);
ALTER table segmentation_feature ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table segmentation_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);


-- feature_type
ALTER table feature_type ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- regulatory_feature
ALTER table regulatory_feature ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);
ALTER table regulatory_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER table regulatory_feature ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);

-- regulatory_attribute
ALTER table regulatory_attribute ADD FOREIGN KEY (regulatory_feature_id) REFERENCES regulatory_feature(regulatory_feature_id);
-- add complex foreign keys

-- annotated_feature
ALTER table annotated_feature ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);
ALTER table annotated_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

-- motif_feature
ALTER table motif_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER table motif_feature ADD FOREIGN KEY (binding_matrix_id) REFERENCES binding_matrix(binding_matrix_id);

-- associated_motif_feature
ALTER table associated_motif_feature ADD FOREIGN KEY (motif_feature_id) REFERENCES motif_feature(motif_feature_id);
ALTER table associated_motif_feature ADD FOREIGN KEY (annotated_feature_id) REFERENCES annotated_feature(annotated_feature_id);

-- binding_matrix
ALTER table binding_matrix ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table binding_matrix ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
-- patched to smallint for v63

-- external_feature
ALTER table external_feature ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);
ALTER table external_feature ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table external_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

-- result_feature
--ALTER table result_feature ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);
--ALTER table result_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
--ERROR 1506 (HY000) at line 72: Foreign key clause is not yet supported in conjunction with partitioning


-- dbfile_registry
-- complex denormalised

-- probe_feature
ALTER table probe_feature ADD FOREIGN KEY (probe_id) REFERENCES probe(probe_id);
ALTER table probe_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table probe_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);


-- data_set
ALTER table data_set ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);

-- supporting_set
ALTER table supporting_set ADD FOREIGN KEY (data_set_id) REFERENCES data_set(data_set_id);
-- + complex denormalised?

-- feature_set
ALTER table feature_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table feature_set ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);
ALTER table feature_set ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table feature_set ADD FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id);


-- result_set
ALTER table result_set ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);
ALTER table result_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table result_set ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- result_set_input
ALTER table result_set_input ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);
-- + complex denormalised ?

-- input_set
ALTER table input_set ADD FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id);
ALTER table input_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table input_set ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);

-- input_subset
ALTER table input_subset ADD FOREIGN KEY (input_set_id) REFERENCES input_set(input_set_id);


-- array_chip
ALTER table array_chip ADD FOREIGN KEY (array_id) REFERENCES array(array_id);

-- probe
ALTER table probe ADD FOREIGN KEY (array_chip_id) REFERENCES array_chip(array_chip_id);
ALTER table probe ADD FOREIGN KEY (probe_set_id) REFERENCES probe_set(probe_set_id);

-- probe_design
ALTER table probe_design ADD FOREIGN KEY (probe_id) REFERENCES probe(probe_id);
ALTER table probe_design ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER table probe_design ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- experiment
ALTER table experiment ADD FOREIGN KEY (experimental_group_id) REFERENCES experimental_group(experimental_group_id);
ALTER table experiment ADD FOREIGN KEY (mage_xml_id) REFERENCES mage_xml(mage_xml_id);

-- experimental_chip
ALTER table experimental_chip ADD FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id);
ALTER table experimental_chip ADD FOREIGN KEY (array_chip_id) REFERENCES array_chip(array_chip_id);
ALTER table experimental_chip ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table experimental_chip ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);

-- channel
ALTER table channel ADD FOREIGN KEY (experimental_chip_id) REFERENCES experimental_chip(experimental_chip_id);

-- result
ALTER table result ADD FOREIGN KEY (result_set_input_id) REFERENCES result_set_input(result_set_input_id);
ALTER table result ADD FOREIGN KEY (probe_id) REFERENCES probe(probe_id);


-- associated_feature_type
ALTER table associated_feature_type ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
-- + complex denormalised


-- Dropping design_type

-- Dropping experimental_design
-- + complex denormalised

-- status

ALTER table status ADD FOREIGN KEY (status_name_id) REFERENCES status_name(status_name_id);
-- + complex denormalised

-- cell_type_lineage
ALTER table cell_type_lineage ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);
ALTER table cell_type_lineage ADD FOREIGN KEY (lineage_id) REFERENCES lineage(lineage_id);



### Core FKs

ALTER table analysis_description ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table external_synonym ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
ALTER table ontology_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER table ontology_xref ADD FOREIGN KEY (source_xref_id) REFERENCES xref(xref_id);
ALTER table object_xref ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table identity_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER table meta_coord ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER table object_xref ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
ALTER table seq_region ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER table unmapped_object ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table unmapped_object ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);
ALTER table unmapped_object ADD FOREIGN KEY (unmapped_reason_id) REFERENCES unmapped_reason(unmapped_reason_id);
ALTER table xref ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);
