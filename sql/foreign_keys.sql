# Foreign key relationships in the Ensembl schema (see table.sql)
#
# This file is intended as a reference since some of the relationships
# are not obvious.
#
# Note that these constraints are not actually used by Ensembl for 
# performance reasons, and referential integrity is enforced at the
# application level. Also MySQL currently does not support foreign
# key constraints on MyISAM tables.


# Core FKs

ALTER table analysis_description ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table external_synonym ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
ALTER table go_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER table go_xref ADD FOREIGN KEY (source_xref_id) REFERENCES xref(xref_id);
ALTER table object_xref ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table identity_xref ADD FOREIGN KEY (object_xref_id) REFERENCES object_xref(object_xref_id);
ALTER table meta_coord ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER table object_xref ADD FOREIGN KEY (xref_id) REFERENCES xref(xref_id);
ALTER table seq_region ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER table unmapped_object ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table unmapped_object ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);
ALTER table unmapped_object ADD FOREIGN KEY (unmapped_reason_id) REFERENCES unmapped_reason(unmapped_reason_id);
ALTER table xref ADD FOREIGN KEY (external_db_id) REFERENCES external_db(external_db_id);


# eFG FKs
ALTER table annotated_feature ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);
ALTER table annotated_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table array_chip ADD FOREIGN KEY (array_id) REFERENCES array(array_id);

ALTER table associated_feature_type ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
-- Not, does not reflect denormalised table_name/id key

ALTER table channel ADD FOREIGN KEY (experimental_chip_id) REFERENCES experimental_chip(experimental_chip_id);

ALTER table result_set_input ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);

ALTER table data_set ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);

-- Dropping design_type

ALTER table experiment ADD FOREIGN KEY (experimental_group_id) REFERENCES experimental_group(experimental_group_id);
ALTER table experiment ADD FOREIGN KEY (mage_xml_id) REFERENCES mage_xml(mage_xml_id);

ALTER table experimental_chip ADD FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id);
ALTER table experimental_chip ADD FOREIGN KEY (array_chip_id) REFERENCES array_chip(array_chip_id);
ALTER table experimental_chip ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table experimental_chip ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);

-- Dropping experimental_design

ALTER table input_set ADD FOREIGN KEY (experiment_id) REFERENCES experiment(experiment_id);
ALTER table input_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table input_set ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);

ALTER table input_subset ADD FOREIGN KEY (input_set_id) REFERENCES input_set(input_set_id);

-- Dropping experimental_variable

ALTER table external_feature ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);
ALTER table external_feature ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table external_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table feature_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table feature_set ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);
ALTER table feature_set ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

ALTER table probe ADD FOREIGN KEY (array_chip_id) REFERENCES array_chip(array_chip_id);
ALTER table probe ADD FOREIGN KEY (probe_set_id) REFERENCES probe_set(probe_set_id);

ALTER table probe_design ADD FOREIGN KEY (probe_id) REFERENCES probe(probe_id);
ALTER table probe_design ADD FOREIGN KEY (coord_system_id) REFERENCES coord_system(coord_system_id);
ALTER table probe_design ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);

ALTER table probe_feature ADD FOREIGN KEY (probe_id) REFERENCES probe(probe_id);
ALTER table probe_feature ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);
ALTER table probe_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table regulatory_attribute ADD FOREIGN KEY (regulatory_feature_id) REFERENCES regulatory_feature(regulatory_feature_id);
-- Can we add complex foreign keys here?

ALTER table regulatory_feature ADD FOREIGN KEY (feature_set_id) REFERENCES feature_set(feature_set_id);
ALTER table regulatory_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);
ALTER table regulatory_feature ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);

ALTER table result ADD FOREIGN KEY (chip_channel_id) REFERENCES chip_channel(chip_channel_id);
ALTER table result ADD FOREIGN KEY (probe_id) REFERENCES probe(probe_id);

ALTER table result_feature ADD FOREIGN KEY (result_set_id) REFERENCES result_set(result_set_id);
ALTER table result_feature ADD FOREIGN KEY (seq_region_id) REFERENCES seq_region(seq_region_id);

ALTER table result_set ADD FOREIGN KEY (cell_type_id) REFERENCES cell_type(cell_type_id);
ALTER table result_set ADD FOREIGN KEY (feature_type_id) REFERENCES feature_type(feature_type_id);
ALTER table result_set ADD FOREIGN KEY (analysis_id) REFERENCES analysis(analysis_id);


ALTER table status ADD FOREIGN KEY (status_name_id) REFERENCES status_name(status_name_id);
--Can we add complex foriegn keys for status here?

ALTER table supporting_set ADD FOREIGN KEY (data_set_id) REFERENCES data_set(data_set_id);
-- Can we had complex foreign keys for different feature tables here?
