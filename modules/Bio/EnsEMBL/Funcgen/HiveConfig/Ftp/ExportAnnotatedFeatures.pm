package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportAnnotatedFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
              '1->A' => 'mk_temp_dir',
              'A->1' => 'rm_annotated_features_temp_dir',
            },
        },
        {   -logic_name  => 'mk_temp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => 'mkdir -p #temp_dir#',
            },
            -flow_into   => {
               MAIN => 'create_chromosome_length_file_for_bigbed_conversion',
            },
        },
        {   -logic_name  => 'create_chromosome_length_file_for_bigbed_conversion',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => 'create_chromosome_length_file_for_bigbed_conversion.pl --registry #reg_conf# --species #species# > #temp_dir#/#species#.ucsc.sizes',
            },
            -flow_into   => {
               MAIN => 'job_factory_annotated_features',
            },
        },
        {   -logic_name  => 'job_factory_annotated_features',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryAnnotatedFeatures',
            -flow_into   => {
               '2->A' => { 'export_annotated_features' => INPUT_PLUS() },
               'A->1' => { 'merge_features'            => INPUT_PLUS() }
            },
        },
        {   -logic_name  => 'export_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 50,
            -batch_size => 100,
            -parameters  => {
	      cmd => 'export_annotated_features.pl --gff_file #temporary_directory#/#partial_gff_file# --bed_file #temporary_directory#/#partial_bed_file# --registry #reg_conf# --species #species# --min_id #min_id# --max_id #max_id# --feature_set_id #feature_set_id#',
            },
        },
        {   -logic_name  => 'merge_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => [
                  'merge_annotated_features_gff',
                  'merge_annotated_features_bed',
                ],
            },
        },
        {   -logic_name  => 'merge_annotated_features_gff',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'cat #gff_files_from_batches# | xargs --max-args=50 cat >> #merged_gff#.unsorted',
            },
            -flow_into   => {
               MAIN => [
                'sort_annotated_features_gff',
                'rm_temp_files_gff',
               ]
            },
        },
        {   -logic_name  => 'rm_temp_files_gff',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'cat #gff_files_from_batches# | xargs --max-args=50 rm -f',
            },
        },
        {   -logic_name  => 'sort_annotated_features_gff',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'sort -k1,1 -k4,4n -o #merged_gff# #merged_gff#.unsorted',
            },
            -flow_into   => {
               MAIN => 'gzip_annotated_features_gff',
            },
        },
        {   -logic_name  => 'gzip_annotated_features_gff',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #merged_gff#',
            },
            -flow_into   => {
               MAIN => 'mk_annotated_features_gff_to_ftp_dir',
            },
        },
        {   -logic_name  => 'mk_annotated_features_gff_to_ftp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mkdir -p #ftp_base_dir#/#species#/AnnotatedFeatures/#epigenome_production_name#/#feature_type_name#',
            },
            -flow_into   => {
               MAIN => 'mv_annotated_features_gff_to_ftp',
            },
        },
        {   -logic_name  => 'mv_annotated_features_gff_to_ftp',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #merged_gff#.gz #ftp_base_dir#/#species#/AnnotatedFeatures/#epigenome_production_name#/#feature_type_name#',
            },
        },
        {   -logic_name  => 'merge_annotated_features_bed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'cat #bed_files_from_batches# | xargs --max-args=50 cat >> #merged_bed#.unsorted',
            },
            -flow_into   => {
               MAIN => [
                'sort_annotated_features_bed',
                'rm_temp_files_bed',
               ]
            },
        },
        {   -logic_name  => 'rm_temp_files_bed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'cat #bed_files_from_batches# | xargs --max-args=50 rm -f',
            },
        },
        {   -logic_name  => 'sort_annotated_features_bed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'sort -k1,1 -k2,2n -o #merged_bed# #merged_bed#.unsorted',
            },
            -flow_into   => {
               MAIN => 'convert_to_bigbed',
            },
        },
        {   -logic_name  => 'convert_to_bigbed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'bedToBigBed #merged_bed# #temp_dir#/#species#.ucsc.sizes #converted_big_bed#',
            },
            -flow_into   => {
               MAIN => 'gzip_annotated_features_bed',
            },
        },
        {   -logic_name  => 'gzip_annotated_features_bed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #merged_bed#',
            },
            -flow_into   => {
               MAIN => 'mk_annotated_features_bed_to_ftp_dir',
            },
        },
        {   -logic_name  => 'mk_annotated_features_bed_to_ftp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mkdir -p #ftp_base_dir#/#species#/AnnotatedFeatures/#epigenome_production_name#/#feature_type_name#',
            },
            -flow_into   => {
               MAIN => [
                'mv_annotated_features_bed_to_ftp',
                'mv_annotated_feature_bigbed_to_ftp'
                ],
            },
        },
        {   -logic_name  => 'mv_annotated_features_bed_to_ftp',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #merged_bed#.gz #ftp_base_dir#/#species#/AnnotatedFeatures/#epigenome_production_name#/#feature_type_name#',
            },
        },
        {   -logic_name  => 'mv_annotated_feature_bigbed_to_ftp',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #converted_big_bed# #ftp_base_dir#/#species#/AnnotatedFeatures/#epigenome_production_name#/#feature_type_name#',
            },
        },
        {   -logic_name  => 'rm_annotated_features_temp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'rm -rf #temp_dir#',
            },
        },
    ]
}

1;
