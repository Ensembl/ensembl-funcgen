-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016] EMBL-European Bioinformatics Institute
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

--  Foreign key relationships in the Ensembl schema (see efg.sql)
--
--  This file is intended as a reference since some of the relationships
--  are not obvious.
--
--  Note that these constraints are not actually used by Ensembl for
--  performance reasons, and referential integrity is enforced at the
--  application level. Also MySQL currently does not support foreign
--  key constraints on MyISAM tables.


-- To run this script you must first create an empty DB and update all tables to InnoDB
-- mysqlw -pXXXXXXX -e 'DROP DATABASE innodb_funcgen'
-- mysqlw -pXXXXXXX -e 'CREATE DATABASE innodb_funcgen'
-- mysqlw -pXXXXXXX innodb_funcgen < efg.sql
-- mysqlro innodb_funcgen -N -e "show tables" | while read t; do if [[ $t != meta ]]; then echo "Altering $t to InnoDB"; mysqlw -pXXXXXXX innodb_funcgen -N -e "ALTER TABLE $t engine=InnoDB"; fi; done
-- Then, optionally disable some to remove complexity when generating ER diagrams, before running:
-- mysqlw -pXXXXXXX innodb_funcgen < foreign_keys.sql

-- Some potential errors
-- ERROR 1071 (42000) at line 1: Specified key was too long; max key length is 767 bytes
-- This is caused by meta TABLE, which has no foreign keys anyway, so skip this.
-- ALTER TABLE   ADD FOREIGN KEY () REFERENCES  ();



-- Last updated for e78_79
ALTER TABLE mirna_target_feature  ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set (feature_set_id);
ALTER TABLE mirna_target_feature  ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type (feature_type_id);


-- Last updated for v75 (included 73,74)
ALTER TABLE result_set  ADD FOREIGN KEY (experiment_id) REFERENCES experiment (experiment_id);
ALTER TABLE feature_set ADD FOREIGN KEY (experiment_id) REFERENCES experiment (experiment_id);


-- 74_75f
ALTER TABLE experiment ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER TABLE experiment ADD FOREIGN KEY (cell_type_id)    REFERENCES cell_type(cell_type_id);

-- 74_75c
ALTER TABLE input_subset ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- 73_74e
-- removed probe_design

-- 73_74_b
-- input_subset
ALTER TABLE input_subset ADD FOREIGN KEY (cell_type_id)    REFERENCES cell_type(cell_type_id);
ALTER TABLE input_subset ADD FOREIGN KEY (experiment_id)   REFERENCES experiment(experiment_id);
ALTER TABLE input_subset ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
-- removed input_subset->input_set link

-- Last updated for v72
ALTER TABLE associated_xref ADD FOREIGN KEY (xref_id)             REFERENCES xref(xref_id);
ALTER TABLE associated_xref ADD FOREIGN KEY (associated_group_id) REFERENCES associated_group(associated_group_id);
ALTER TABLE associated_xref ADD FOREIGN KEY (object_xref_id)      REFERENCES identity_xref(object_xref_id);

-- Removed 75_76_b
-- feature_set
-- ALTER TABLE feature_set ADD FOREIGN KEY (input_set_id) REFERENCES input_set(input_set_id);

-- input_subset_tracking
-- ALTER TABLE input_subset_tracking ADD FOREIGN KEY (input_subset_id) REFERENCES input_subset(input_subset_id);


-- input_set_tracking
-- ALTER TABLE input_set_tracking ADD FOREIGN KEY (input_set_id) REFERENCES input_set(input_set_id);


-- result_set_tracking
-- ALTER TABLE result_set_tracking ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);


-- data_set_tracking
-- ALTER TABLE data_set_tracking ADD FOREIGN KEY (data_set_id) REFERENCES data_set(data_set_id);

### eFG FKs

-- segmentation_feature
ALTER TABLE segmentation_feature ADD FOREIGN KEY (feature_set_id)   REFERENCES feature_set(feature_set_id);
ALTER TABLE segmentation_feature ADD FOREIGN KEY (feature_type_id)  REFERENCES feature_type(feature_type_id);
ALTER TABLE segmentation_feature ADD FOREIGN KEY (seq_region_id)    REFERENCES seq_region(seq_region_id);


-- feature_type
ALTER TABLE feature_type ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- regulatory_feature
ALTER TABLE regulatory_feature ADD FOREIGN KEY (feature_set_id)   REFERENCES feature_set(feature_set_id);
ALTER TABLE regulatory_feature ADD FOREIGN KEY (seq_region_id)    REFERENCES seq_region(seq_region_id);
ALTER TABLE regulatory_feature ADD FOREIGN KEY (feature_type_id)  REFERENCES feature_type(feature_type_id);

-- regulatory_attribute
ALTER TABLE regulatory_attribute ADD FOREIGN KEY (regulatory_feature_id) REFERENCES regulatory_feature(regulatory_feature_id);
-- add complex foreign keys?

-- annotated_feature
ALTER TABLE annotated_feature ADD FOREIGN KEY (feature_set_id)  REFERENCES feature_set(feature_set_id);
ALTER TABLE annotated_feature ADD FOREIGN KEY (seq_region_id)   REFERENCES seq_region(seq_region_id);

-- motif_feature
ALTER TABLE motif_feature ADD FOREIGN KEY (seq_region_id)     REFERENCES seq_region(seq_region_id);
ALTER TABLE motif_feature ADD FOREIGN KEY (binding_matrix_id) REFERENCES binding_matrix(binding_matrix_id);

-- associated_motif_feature
ALTER TABLE associated_motif_feature ADD FOREIGN KEY (motif_feature_id)     REFERENCES motif_feature(motif_feature_id);
ALTER TABLE associated_motif_feature ADD FOREIGN KEY (annotated_feature_id) REFERENCES annotated_feature(annotated_feature_id);

-- binding_matrix
ALTER TABLE binding_matrix ADD FOREIGN KEY (feature_type_id)  REFERENCES feature_type(feature_type_id);
ALTER TABLE binding_matrix ADD FOREIGN KEY (analysis_id)      REFERENCES analysis(analysis_id);
-- patched to smallint for v63

-- external_feature
ALTER TABLE external_feature ADD FOREIGN KEY (feature_set_id)   REFERENCES feature_set(feature_set_id);
ALTER TABLE external_feature ADD FOREIGN KEY (feature_type_id)  REFERENCES feature_type(feature_type_id);
ALTER TABLE external_feature ADD FOREIGN KEY (seq_region_id)    REFERENCES seq_region(seq_region_id);

-- result_feature
ALTER TABLE result_feature ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);
ALTER TABLE result_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);



-- dbfile_registry
-- complex denormalised

-- probe_feature
ALTER TABLE probe_feature ADD FOREIGN KEY (probe_id)      REFERENCES probe(probe_id);
ALTER TABLE probe_feature ADD FOREIGN KEY (analysis_id)   REFERENCES analysis(analysis_id);
ALTER TABLE probe_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);


-- data_set
ALTER TABLE data_set ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);

-- supporting_set
ALTER TABLE supporting_set ADD FOREIGN KEY (data_set_id) REFERENCES data_set(data_set_id);
-- + complex denormalised?

-- feature_set
ALTER TABLE feature_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER TABLE feature_set ADD FOREIGN KEY (cell_type_id)    REFERENCES cell_type(cell_type_id);
ALTER TABLE feature_set ADD FOREIGN KEY (analysis_id)     REFERENCES analysis(analysis_id);


-- result_set
ALTER TABLE result_set ADD FOREIGN KEY (cell_type_id)    REFERENCES cell_type(cell_type_id);
ALTER TABLE result_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER TABLE result_set ADD FOREIGN KEY (analysis_id)     REFERENCES analysis(analysis_id);

-- result_set_input
ALTER TABLE result_set_input ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);
-- + complex denormalised ?

-- input_set
ALTER TABLE input_set ADD FOREIGN KEY (experiment_id)   REFERENCES experiment(experiment_id);
ALTER TABLE input_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER TABLE input_set ADD FOREIGN KEY (cell_type_id)    REFERENCES cell_type(cell_type_id);

-- input_subset
-- removed r75
-- ALTER TABLE input_subset ADD FOREIGN KEY (input_set_id) REFERENCES input_set(input_set_id);


-- array_chip
ALTER TABLE array_chip ADD FOREIGN KEY (array_id) REFERENCES array(array_id);

-- probe
ALTER TABLE probe ADD FOREIGN KEY (array_chip_id) REFERENCES array_chip(array_chip_id);
ALTER TABLE probe ADD FOREIGN KEY (probe_set_id)  REFERENCES probe_set(probe_set_id);

-- probe_design
-- removed 73_74e
-- ALTER TABLE probe_design ADD FOREIGN KEY (probe_id) REFERENCES probe(probe_id);
-- ALTER TABLE probe_design ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
-- ALTER TABLE probe_design ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- experiment
ALTER TABLE experiment ADD FOREIGN KEY (experimental_group_id)  REFERENCES experimental_group(experimental_group_id);
ALTER TABLE experiment ADD FOREIGN KEY (mage_xml_id)            REFERENCES mage_xml(mage_xml_id);

-- experimental_chip
ALTER TABLE experimental_chip ADD FOREIGN KEY (experiment_id)   REFERENCES experiment(experiment_id);
ALTER TABLE experimental_chip ADD FOREIGN KEY (array_chip_id)   REFERENCES array_chip(array_chip_id);
ALTER TABLE experimental_chip ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER TABLE experimental_chip ADD FOREIGN KEY (cell_type_id)    REFERENCES cell_type(cell_type_id);

-- channel
ALTER TABLE channel ADD FOREIGN KEY (experimental_chip_id) REFERENCES experimental_chip(experimental_chip_id);

-- result
ALTER TABLE result ADD FOREIGN KEY (result_set_input_id)  REFERENCES result_set_input(result_set_input_id);
ALTER TABLE result ADD FOREIGN KEY (probe_id)             REFERENCES probe(probe_id);


-- associated_feature_type
ALTER TABLE associated_feature_type ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
-- + complex denormalised


-- Dropping design_type

-- Dropping experimental_design
-- + complex denormalised

-- status

ALTER TABLE status ADD FOREIGN KEY (status_name_id) REFERENCES status_name(status_name_id);
-- + complex denormalised

-- cell_type_lineage
ALTER TABLE cell_type_lineage ADD FOREIGN KEY (cell_type_id)  REFERENCES cell_type(cell_type_id);
ALTER TABLE cell_type_lineage ADD FOREIGN KEY (lineage_id)    REFERENCES lineage(lineage_id);



-- regbuild_string
-- None



### Core FKs

ALTER TABLE analysis_description ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE external_synonym ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
ALTER TABLE ontology_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER TABLE ontology_xref ADD FOREIGN KEY (source_xref_id) REFERENCES xref(xref_id);
ALTER TABLE object_xref ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE identity_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER TABLE meta_coord ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER TABLE object_xref ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
ALTER TABLE seq_region ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER TABLE unmapped_object ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE unmapped_object ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);
ALTER TABLE unmapped_object ADD FOREIGN KEY (unmapped_reason_id) REFERENCES unmapped_reason(unmapped_reason_id);
ALTER TABLE xref ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);
