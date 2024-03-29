Ensembl Funcgen Deprecated Methods
==================================

This file contains the list of methods, modules and scripts deprecated in the Ensembl Funcgen API.
A method is deprecated when it is not functional any more (schema/data change) or has been replaced by a better one.
When a method is deprecated, a deprecation warning is thrown whenever the method is used.

### To be removed in Ensembl Release 109 ###
- Bio::EnsEMBL::Funcgen::**MotifFeature**::*get_all_overlapping_Peaks()*
- Bio::EnsEMBL::Funcgen::**MotifFeature**::*get_all_overlapping_Peaks_by_Epigenome()*
- Bio::EnsEMBL::Funcgen::DBSQL::**MotifFeatureAdaptor**::*_fetch_all_overlapping_Peaks()*
- Bio::EnsEMBL::Funcgen::DBSQL::**MotifFeatureAdaptor**::*_fetch_all_overlapping_Peaks_by_Epigenome()*
- Bio::EnsEMBL::Funcgen::DBSQL::**MotifFeatureAdaptor**::*store_associated_Peak()*
- Bio::EnsEMBL::Funcgen::**Peak**::*_constructor_parameters()*
- Bio::EnsEMBL::Funcgen::**Peak**::*dbID()*
- Bio::EnsEMBL::Funcgen::**Peak**::*adaptor()*
- Bio::EnsEMBL::Funcgen::**Peak**::*peak_calling_id()*
- Bio::EnsEMBL::Funcgen::**Peak**::*summit()*
- Bio::EnsEMBL::Funcgen::**Peak**::*score()*
- Bio::EnsEMBL::Funcgen::**Peak**::*start()*
- Bio::EnsEMBL::Funcgen::**Peak**::*end()*
- Bio::EnsEMBL::Funcgen::**Peak**::*seq_region_id()*
- Bio::EnsEMBL::Funcgen::**Peak**::*seq_region_start()*
- Bio::EnsEMBL::Funcgen::**Peak**::*seq_region_end()*
- Bio::EnsEMBL::Funcgen::**Peak**::*seq_region_strand()*
- Bio::EnsEMBL::Funcgen::**Peak**::*strand()*
- Bio::EnsEMBL::Funcgen::**Peak**::*slice()*
- Bio::EnsEMBL::Funcgen::**Peak**::*get_PeakCalling()*
- Bio::EnsEMBL::Funcgen::**Peak**::*set_PeakCalling()*
- Bio::EnsEMBL::Funcgen::**Peak**::*display_label()*
- Bio::EnsEMBL::Funcgen::**Peak**::*display_id()*
- Bio::EnsEMBL::Funcgen::**Peak**::*get_all_MotifFeatures()*
- Bio::EnsEMBL::Funcgen::**Peak**::*seq_region_name()*
- Bio::EnsEMBL::Funcgen::**Peak**::*feature_so_acc()*
- Bio::EnsEMBL::Funcgen::**Peak**::*feature_so_term()*
- Bio::EnsEMBL::Funcgen::**Peak**::*summary_as_hash()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*object_class()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*_tables()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*insertion_method()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*_columns()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*fetch_all_by_PeakCalling()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*fetch_all_by_Slice_PeakCalling()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*_fetch_overlapping_MotifFeatures()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*_parse_bed_line()*
- Bio::EnsEMBL::Funcgen::DBSQL::**PeakAdaptor**::*_bulk_export_to_bed_by_PeakCalling()*


### Removed in EnsEMBL Release 104 ###
 - Bio::EnsEMBL::Funcgen::**FeatureType**::*so_name()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_Analysis()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_Chance()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_control_Alignment()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_Epigenome()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_Experiment()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_FeatureType()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_Frip()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_Idr()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_PeakCallingStatistic()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_signal_Alignment()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*fetch_source_label()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_all_deduplicated_replicate_Alignments()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_all_ReadFileExperimentalConfigurations()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_all_ReadFiles()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_Analysis()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_bam_DataFile()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_bigwig_DataFile()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_Chance_by_control_Alignment()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_Experiment()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_PhantomPeak()*
 - Bio::EnsEMBL::Funcgen::**Alignment**::*fetch_source_Alignment()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryBuild**:*fetch_Analysis()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryBuild**:*fetch_FeatureType()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryBuild**:*fetch_sample_RegulatoryFeature()*
 - Bio::EnsEMBL::Funcgen::**SegmentationStateEmission**:*fetch_segmentation_state_assignment()*
 - Bio::EnsEMBL::Funcgen::**Frip**:*fetch_PeakCalling()*
 - Bio::EnsEMBL::Funcgen::**Idr**:*fetch_Experiment()*
 - Bio::EnsEMBL::Funcgen::**ProbeTranscriptMapping**:*fetch_Probe()*
 - Bio::EnsEMBL::Funcgen::**ProbeSetTranscriptMapping**:*fetch_ProbeSet()*
 - Bio::EnsEMBL::Funcgen::**ProbeFeatureTranscriptMapping**:*fetch_ProbeFeature()*
 - Bio::EnsEMBL::Funcgen::**PeakCallingStatistic**:*fetch_PeakCalling()*
 - Bio::EnsEMBL::Funcgen::**PhantomPeak**:*fetch_Alignment()*
 - Bio::EnsEMBL::Funcgen::**PhantomPeak**:*fetch_Idr()*
 - Bio::EnsEMBL::Funcgen::**ReadFile**::*fetch_FastQC()*
 - Bio::EnsEMBL::Funcgen::**ReadFile**::*fetch_mate_ReadFile()*
 - Bio::EnsEMBL::Funcgen::**ReadFile**::*fetch_ReadFileExperimentalConfiguration()*
 - Bio::EnsEMBL::Funcgen::**Probe**::*fetch_all_ProbeTranscriptMappings()*
 - Bio::EnsEMBL::Funcgen::**ProbeSet**:*fetch_all_ProbeSetTranscriptMappings()*
 - Bio::EnsEMBL::Funcgen::**Chance**::*db()*
 - Bio::EnsEMBL::Funcgen::**DataFile**::*db()*
 - Bio::EnsEMBL::Funcgen::**ExampleFeature**::*db()*
 - Bio::EnsEMBL::Funcgen::**ExecutionPlan**::*db()*
 - Bio::EnsEMBL::Funcgen::**FastQC**::*db()*
 - Bio::EnsEMBL::Funcgen::**Frip**::*db()*
 - Bio::EnsEMBL::Funcgen::**Idr**::*db()*
 - Bio::EnsEMBL::Funcgen::**Peak**::*db()*
 - Bio::EnsEMBL::Funcgen::**PeakCalling**::*db()*
 - Bio::EnsEMBL::Funcgen::**PeakCallingStatistic**::*db()*
 - Bio::EnsEMBL::Funcgen::**PhantomPeak**::*db()*
 - Bio::EnsEMBL::Funcgen::**ProbeMapping**::*db()*
 - Bio::EnsEMBL::Funcgen::**ProbeMappingStatistic**::*db()*
 - Bio::EnsEMBL::Funcgen::**ReadFile**::*db()*
 - Bio::EnsEMBL::Funcgen::**ReadFileExperimentalConfiguration**::*db()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryActivity**::*db()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryBuildStatistic**::*db()*
 - Bio::EnsEMBL::Funcgen::**Segmentation**::*db()*
 - Bio::EnsEMBL::Funcgen::**SegmentationStateAssignment**::*db()*
 - Bio::EnsEMBL::Funcgen::**SegmentationStateEmission**::*db()*
 - Bio::EnsEMBL::Funcgen::**SegmentationStatistic**::*db()*

### Removed in EnsEMBL Release 101 ###

 - Bio::EnsEMBL::Funcgen::**BindingMatrixFrequencies**::*binding_matrix*
 - Bio::EnsEMBL::Funcgen::**MotifFeature**::*binding_matrix*
 - Bio::EnsEMBL::Funcgen::DBSQL::**EpigenomeAdaptor**::*fetch_by_display_label()*
 - Bio::EnsEMBL::Funcgen::**Epigenome**::*display_label()*

### Removed in EnsEMBL Release 100 ###

 - Bio::EnsEMBL::Funcgen::**MotifFeature**::*fetch_overlapping_Peak_by_Epigenome()*
 - Bio::EnsEMBL::Funcgen::**MotifFeature**::*fetch_all_overlapping_Peaks()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**MotifFeatureAdaptor**::*fetch_all_by_Slice_Epigenome()*
 - Bio::EnsEMBL::Funcgen::**Peak**::*fetch_PeakCalling()*
 - Bio::EnsEMBL::Funcgen::**Peak**::*get_underlying_structure()*
 - Bio::EnsEMBL::Funcgen::**Peak**::*fetch_all_MotifFeatures()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*fetch_all_MotifFeatures()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*fetch_all_MotifFeatures_with_matching_Peak()*

### Removed in EnsEMBL Release 99 ###

 - Bio::EnsEMBL::Funcgen::**Alignment**::*db*

### Removed in EnsEMBL Release 97 ###

 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*fetch_all_MotifFeatures_by_Epigenome()*

### Removed in EnsEMBL Release 94 ###

 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*feature_type()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*description()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*analysis()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*frequencies()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*_build_matrix()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*_validate_matrix()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*frequency_matrix()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*weight_matrix()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*matrix()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*weights()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*_process_frequency_matrix()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*_compute_base_weight()*
 - Bio::EnsEMBL::Funcgen::**BindingMatrix**::*info_content()*
 - Bio::EnsEMBL::Funcgen::**MotifFeature**::*feature_type()*
 - Bio::EnsEMBL::Funcgen::**MotifFeature**::*display_label()*
 - Bio::EnsEMBL::Funcgen::**MotifFeature**::*associated_annotated_features()*
 - Bio::EnsEMBL::Funcgen::**MotifFeature**::*interdb_stable_id()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**BindingMatrixAdaptor**::*fetch_all_by_name()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**BindingMatrixAdaptor**::*fetch_all_by_name_FeatureType()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**BindingMatrixAdaptor**::*fetch_all_by_FeatureType()*

### Removed in Ensembl Release 93 ###

 - Bio::EnsEMBL::Funcgen::Epigenome::ontology_accession()
 - Bio::EnsEMBL::Funcgen::Epigenome::tissue()
 - Bio::EnsEMBL::Funcgen::Probe::get_complete_name()
 - Bio::EnsEMBL::Funcgen::RegulatoryFeature::is_projected()

### Removed in Ensembl Release 92 ###

 - Bio::EnsEMBL::Funcgen::**RegulatoryActivity**::*regulatory_evidence()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*regulatory_evidence()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**ProbeAdaptor**::*fetch_all_by_external_name()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**ProbeSetAdaptor**::*fetch_all_by_external_name()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**ProbeSetAdaptor**::*fetch_by_array_probeset_name()*
 - Bio::EnsEMBL::Funcgen::**Probe**::*get_all_Transcript_DBEntries()*
 - Bio::EnsEMBL::Funcgen::**ProbeSetTranscriptMapping**::*display_id()*
 - Bio::EnsEMBL::Funcgen::**ProbeSetTranscriptMapping**::*linkage_annotation()*
 - Bio::EnsEMBL::Funcgen::**ProbeTranscriptMapping**::*display_id()*
 - Bio::EnsEMBL::Funcgen::**ProbeTranscriptMapping**::*linkage_annotation()*

### Removed in Ensembl Release 90 ###
 - Bio::EnsEMBL::Funcgen::DBSQL::**RegulatoryFeatureAdaptor**::*fetch_all_by_Slice_Epigenomes()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**RegulatoryFeatureAdaptor**::*fetch_all_by_Slice_Epigenomes_Activity()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**RegulatoryFeatureAdaptor**::*fetch_all_by_Slice_Activity()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**FeatureSetAdaptor**::*fetch_attribute_set_config_by_FeatureSet()*
 - Bio::EnsEMBL::Funcgen::**FeatureSet**::*is_attribute_set()*
 - Bio::EnsEMBL::Funcgen::Parsers::**InputSet**
 - Bio::EnsEMBL::Funcgen::Parsers::**Bed**
 - Bio::EnsEMBL::Funcgen::Parsers::**GFF**
 - Bio::EnsEMBL::Funcgen::Parsers::**InputSet**
 - Bio::EnsEMBL::Funcgen::Parsers::**MAGE**
 - Bio::EnsEMBL::Funcgen::Parsers::**Nimblegen** 

### Removed in Ensembl Release 89 ###
 - Bio::EnsEMBL::Funcgen::**ResultSet**::*replicate()*
 - Bio::Ensembl::Funcgen::DBSQL::**CellTypeAdaptor**
 - Bio::Ensembl::Funcgen::**CellType**
 - Bio::EnsEMBL::Funcgen::**DataSet**::*cell_type()*
 - Bio::EnsEMBL::Funcgen::**Experiment**::*cell_type()*
 - Bio::EnsEMBL::Funcgen::**Importer**::*cell_type()*
 - Bio::EnsEMBL::Funcgen::**SetFeature**::*cell_type()*
 - Bio::EnsEMBL::Funcgen::**Set**::*cell_type()*
 - Bio::EnsEMBL::Funcgen::**Epigenome**::*efo_id()*
 - Bio::EnsEMBL::Funcgen::**Experiment**::*primary_design_type()*
 - Bio::EnsEMBL::Funcgen::**Experiment**::*description()*
 - Bio::EnsEMBL::Funcgen::**Experiment**::*mage_xml()*
 - Bio::EnsEMBL::Funcgen::**Experiment**::*mage_xml_id()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*cell_type_count()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**ExperimentAdaptor**::*fetch_all_by_CellType()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**MotifFeatureAdaptor**::*fetch_all_by_Slice_CellType()*
 - Bio::EnsEMBL::Funcgen::DBSQL::**SetAdaptor**::*fetch_all_by_CellType()*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*get_focus_attributes*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*get_nonfocus_attributes*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*is_unique_to_FeatureSets*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeature**::*get_other_RegulatoryFeatures*
 - Bio::EnsEMBL::Funcgen::**RegulatoryFeatureAdaptor**::*fetch_all_by_stable_ID*

### Removed in Ensembl Release 88 ###
 - Bio::Ensembl::Funcgen::**Channel**
 - Bio::Ensembl::Funcgen::DBSQL::**ChannelAdaptor**
 - Bio::Ensembl::Funcgen::**ExperimentalChip**
 - Bio::Ensembl::Funcgen::DBSQL::**ExperimentalChipAdaptor**
 - Bio::Ensembl::Funcgen::**InputSet**
 - Bio::Ensembl::Funcgen::DBSQL::**InputSetAdaptor**
 - Bio::Ensembl::Funcgen::**ResultFeature**
 - Bio::Ensembl::Funcgen::DBSQL::**ResultFeatureAdaptor**
 - Bio::Ensembl::Funcgen::Collector::**ResultFeature**

### Removed in Ensembl Release 86 ###
 - Bio::Ensembl::Funcgen::**ResultSet**::*get_ResultFeatures_by_Slice*

### Removed in Ensembl Release 85 ###
 - Bio::Ensembl::Funcgen::**Experiment**::*date()*
 - Bio::Ensembl::Funcgen::**ResultSet**::*get_replicate_set_by_chip_channel_id()*
 - Bio::Ensembl::Funcgen::DBSQL::**BaseAdaptor**::*list_dbIDs()*
 - Bio::Ensembl::Funcgen::DBSQL::**BaseAdaptor**::*_constrain_status()*
 - Bio::Ensembl::Funcgen::DBSQL::**BaseAdaptor**::*fetch_all_by_status()*
 - Bio::Ensembl::Funcgen::DBSQL::**DBEntryAdaptor**::*list_regulatory_feature_ids_by_extid()*
 - Bio::Ensembl::Funcgen::DBSQL::**DBEntryAdaptor**::*list_probeset_ids_by_extid()*
 - Bio::Ensembl::Funcgen::DBSQL::**DBEntryAdaptor**::*list_feature_type_ids_by_extid()*
 - Bio::Ensembl::Funcgen::DBSQL::**DBEntryAdaptor**::*list_external_feature_ids_by_extid()*
 - Bio::Ensembl::Funcgen::DBSQL::**DBEntryAdaptor**::*list_annotated_feature_ids_by_extid()*
 - Bio::Ensembl::Funcgen::DBSQL::**DBEntryAdaptor**::*list_probe_ids_by_extid()*
 - Bio::Ensembl::Funcgen::DBSQL::**FeatureSetAdaptor**::*fetch_all_by_type()*
 - Bio::Ensembl::Funcgen::DBSQL::**InputSubsetAdaptor**::*fetch_by_name_and_experiment()*
 - Bio::Ensembl::Funcgen::DBSQL::**ProbeFeatureAdaptor**::*fetch_all_by_probeset()*
 - Bio::Ensembl::Funcgen::DBSQL::**ResultFeatureAdaptor**::*fetch_all()*
 - Bio::Ensembl::Funcgen::DBSQL::**ResultFeatureAdaptor**::*fetch_by_dbID()*
 - Bio::Ensembl::Funcgen::DBSQL::**ResultFeatureAdaptor**::*fetch_all_by_dbID_list()*
 - Bio::Ensembl::Funcgen::DBSQL::**ResultFeatureAdaptor**::*fetch_all_by_logic_name()*
 - Bio::Ensembl::Funcgen::DBSQL::**ResultFeatureAdaptor**::*_list_seq_region_ids()*
 - Bio::Ensembl::Funcgen::DBSQL::**ResultSetAdaptor**::*store_dbfile_data_dir()*
 - Bio::Ensembl::Funcgen::DBSQL::**ResultSetAdaptor**::*_fetch_dbfile_data_dir()*

### Removed in Ensembl Release 84 ###
 - Bio::Ensembl::Funcgen::**Dataset**::*add_supporting_sets()*
 - Bio::Ensembl::Funcgen::**Dataset**::*_validate_and_set_types()*
 - Bio::Ensembl::Funcgen::**InputSubset**::*input_set()*
 - Bio::Ensembl::Funcgen::**InputSubset**::*archive_id()*
 - Bio::Ensembl::Funcgen::**InputSubset**::*display_url()*
 - Bio::Ensembl::Funcgen::**Probe**::*add_Analysis_score()*
 - Bio::Ensembl::Funcgen::**Probe**::*add_Analysis_CoordSystem_score()*
 - Bio::Ensembl::Funcgen::**Probe**::*get_score_by_Analysis()*
 - Bio::Ensembl::Funcgen::**Probe**::*get_score_by_Analysis_CoordSystem()*
 - Bio::Ensembl::Funcgen::**Probe**::*get_all_design_scores()*
 - Bio::Ensembl::Funcgen::DBSQL::**DataSetAdaptor**::*store_updated_sets()*
 - Bio::Ensembl::Funcgen::HiveConfig::**Alignment_conf.pm**
 - Bio::Ensembl::Funcgen::HiveConfig::**Peaks_conf.pm**
 - Bio::Ensembl::Funcgen::HiveConfig::**Annotation_conf.pm**
 - Bio::Ensembl::Funcgen::HiveConfig::**Dnase_profile_conf.pm**
 - Bio::Ensembl::Funcgen::HiveConfig::**MotifFinder_conf.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**Alignment.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**AnnotateRegulatoryFeatures.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**Annotation.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**ClusterMotifs.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**ConvergeReplicates.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**Funcgen.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**Import.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**InferMotifs.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**MakeDnaseProfile.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**Motif.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**RunBWA.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**RunCCAT.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**RunCentipede.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**RunMACS.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**RunSWEmbl.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**SWEmbl.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**SetupAlignmentPipeline.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**SetupAnnotationPipeline.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**SetupMotifInference.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**SetupPeaksPipeline.pm**
 - Bio::EnsEMBL::Funcgen::RunnableDB::**WrapUpAlignment.pm**

 - scripts/pipeline/**configure_hive.pl**
 - scripts/regulatory_build/**load_segmentation.pl**
 - scripts/regulatory_build/**unpack_swilder_segmentation.pl**
