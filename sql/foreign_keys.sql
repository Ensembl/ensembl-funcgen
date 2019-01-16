-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
-- mysql -pXXXXXXX -e 'DROP DATABASE innodb_funcgen'
-- mysql -pXXXXXXX -e 'CREATE DATABASE innodb_funcgen'
-- mysql -pXXXXXXX innodb_funcgen < table.sql
-- mysql  innodb_funcgen -N -e "show tables" | while read t; do if [[ $t != meta ]]; then echo "Altering $t to InnoDB"; mysqlw -pXXXXXXX innodb_funcgen -N -e "ALTER TABLE $t engine=InnoDB"; fi; done
-- Then, optionally disable some to remove complexity when generating ER diagrams, before running:
-- mysqlw -pXXXXXXX innodb_funcgen < foreign_keys.sql

-- Some potential errors
-- ERROR 1071 (42000) at line 1: Specified key was too long; max key length is 767 bytes
-- This is caused by meta TABLE, which has no foreign keys anyway, so skip this.
-- ALTER TABLE   ADD FOREIGN KEY () REFERENCES  ();


-- regulatory_feature
ALTER TABLE regulatory_feature ADD FOREIGN KEY (feature_type_id)  REFERENCES feature_type(feature_type_id);
ALTER TABLE regulatory_feature ADD FOREIGN KEY (regulatory_build_id)   REFERENCES regulatory_build(regulatory_build_id);

-- regulatory_evidence
ALTER TABLE regulatory_evidence ADD FOREIGN KEY (regulatory_activity_id) REFERENCES regulatory_activity(regulatory_activity_id);
	-- TODO polymorphic association pending

-- regulatory_activity
ALTER TABLE regulatory_activity ADD FOREIGN KEY (regulatory_feature_id) REFERENCES regulatory_feature(regulatory_feature_id);
ALTER TABLE regulatory_activity ADD FOREIGN KEY (epigenome_id) REFERENCES epigenome(epigenome_id);

-- regulatory_build
ALTER TABLE regulatory_build ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER TABLE regulatory_build ADD FOREIGN KEY (analysis_id) REFERENCES analysis (analysis_id);
ALTER TABLE regulatory_build ADD FOREIGN KEY (sample_regulatory_feature_id) REFERENCES regulatory_feature(regulatory_feature_id);

-- regulatory_build_epigenome
ALTER TABLE regulatory_build_epigenome ADD FOREIGN KEY (regulatory_build_id) REFERENCES regulatory_build(regulatory_build_id);
ALTER TABLE regulatory_build_epigenome ADD FOREIGN KEY (epigenome_id) REFERENCES epigenome(epigenome_id);

-- segmentation_file
ALTER TABLE segmentation_file ADD FOREIGN KEY (regulatory_build_id) REFERENCES regulatory_build(regulatory_build_id);
ALTER TABLE segmentation_file ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE segmentation_file ADD FOREIGN KEY (epigenome_id) REFERENCES epigenome(epigenome_id);

-- motif_feature
ALTER TABLE motif_feature ADD FOREIGN KEY (binding_matrix_id) REFERENCES binding_matrix(binding_matrix_id);

-- mirna_target_feature
ALTER TABLE mirna_target_feature  ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set (feature_set_id);
ALTER TABLE mirna_target_feature  ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type (feature_type_id);

-- associated_motif_feature
ALTER TABLE associated_motif_feature ADD FOREIGN KEY (motif_feature_id) REFERENCES motif_feature(motif_feature_id);

-- external_feature
ALTER TABLE external_feature ADD FOREIGN KEY (feature_set_id)   REFERENCES feature_set(feature_set_id);
ALTER TABLE external_feature ADD FOREIGN KEY (feature_type_id)  REFERENCES feature_type(feature_type_id);

-- external_feature_file
ALTER TABLE external_feature_file ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE external_feature_file ADD FOREIGN KEY (epigenome_id) REFERENCES epigenome(epigenome_id);
ALTER TABLE external_feature_file ADD FOREIGN KEY (feature_type_id)  REFERENCES feature_type(feature_type_id);

-- probe_feature
ALTER TABLE probe_feature ADD FOREIGN KEY (probe_id)      REFERENCES probe(probe_id);
ALTER TABLE probe_feature ADD FOREIGN KEY (analysis_id)   REFERENCES analysis(analysis_id);

-- feature_type
ALTER TABLE feature_type ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- associated_feature_type
ALTER TABLE associated_feature_type ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
	-- TODO polymorphic association pending

-- feature_set
ALTER TABLE feature_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER TABLE feature_set ADD FOREIGN KEY (analysis_id)     REFERENCES analysis(analysis_id);

-- dbfile_registry
	-- TODO polymorphic association pending

-- array_chip
ALTER TABLE array_chip ADD FOREIGN KEY (array_id) REFERENCES array(array_id);

-- probe
ALTER TABLE probe ADD FOREIGN KEY (array_chip_id) REFERENCES array_chip(array_chip_id);
ALTER TABLE probe ADD FOREIGN KEY (probe_set_id)  REFERENCES probe_set(probe_set_id);
ALTER TABLE probe ADD FOREIGN KEY (probe_seq_id)  REFERENCES probe_seq(probe_seq_id);

-- probe_feature_transcript
ALTER TABLE probe_feature_transcript ADD FOREIGN KEY (probe_feature_id) REFERENCES probe_feature(probe_feature_id);

-- probe_set
ALTER TABLE probe_set ADD FOREIGN KEY (array_chip_id) REFERENCES array_chip(array_chip_id);

-- probe_set_transcript
ALTER TABLE probe_set_transcript ADD FOREIGN KEY (probe_set_id) REFERENCES probe_set(probe_set_id);

-- experiment
ALTER TABLE experiment ADD FOREIGN KEY (experimental_group_id)  REFERENCES experimental_group(experimental_group_id);
ALTER TABLE experiment ADD FOREIGN KEY (control_id) REFERENCES experiment(experiment_id);
ALTER TABLE experiment ADD FOREIGN KEY (epigenome_id)    REFERENCES epigenome(epigenome_id);
ALTER TABLE experiment ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);

-- alignment
ALTER TABLE alignment ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE alignment ADD FOREIGN KEY (bam_file_id) REFERENCES data_file(data_file_id);
ALTER TABLE alignment ADD FOREIGN KEY (bigwig_file_id) REFERENCES data_file(data_file_id);

-- alignment_read_file
ALTER TABLE alignment_read_file ADD FOREIGN KEY (alignment_id) REFERENCES alignment(alignment_id);
ALTER TABLE alignment_read_file ADD FOREIGN KEY (read_file_id) REFERENCES read_file(read_file_id);

-- peak
ALTER TABLE peak ADD FOREIGN KEY (peak_calling_id) REFERENCES peak_calling(peak_calling_id);

-- read_file
ALTER TABLE read_file ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- read_file_experimental_configuration
ALTER TABLE read_file_experimental_configuration ADD FOREIGN KEY (read_file_id) REFERENCES read_file(read_file_id);
ALTER TABLE read_file_experimental_configuration ADD FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id);

-- analysis_description
ALTER TABLE analysis_description ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

-- associated_xref
ALTER TABLE associated_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER TABLE associated_xref ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
ALTER TABLE associated_xref ADD FOREIGN KEY (source_xref_id) REFERENCES xref(xref_id);
ALTER TABLE associated_xref ADD FOREIGN KEY (associated_group_id) REFERENCES associated_group(associated_group_id);

-- identity_xref
ALTER TABLE identity_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);

-- external_synonym
ALTER TABLE external_synonym ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);

-- ontology_xref
ALTER TABLE ontology_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER TABLE ontology_xref ADD FOREIGN KEY (source_xref_id) REFERENCES xref(xref_id);

-- xref
ALTER TABLE xref ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);

-- object_xref
ALTER TABLE object_xref ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE object_xref ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
	-- TODO polymorphic association pending


-- unmapped_object
ALTER TABLE unmapped_object ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER TABLE unmapped_object ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);
ALTER TABLE unmapped_object ADD FOREIGN KEY (unmapped_reason_id) REFERENCES unmapped_reason(unmapped_reason_id);
	-- TODO polymorphic association pending


-- TRACKING TABLES
-- input_subset_tracking
-- ALTER TABLE input_subset_tracking ADD FOREIGN KEY (input_subset_id) REFERENCES input_subset(input_subset_id);


-- result_set_tracking
-- ALTER TABLE result_set_tracking ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);


-- data_set_tracking
-- ALTER TABLE data_set_tracking ADD FOREIGN KEY (data_set_id) REFERENCES data_set(data_set_id);



