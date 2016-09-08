package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportAnnotatedFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      ftp_base_dir  => $self->o('ftp_base_dir'),
      reg_conf      => $self->o('reg_conf'),
      temp_dir      => '#ftp_base_dir#/#species#/tempdir/annotated_features',
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN' => 'job_factory_annotated_features'
            },
        },
        {   -logic_name  => 'job_factory_annotated_features',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryAnnotatedFeatures',
            -flow_into   => {
               '2->A' => { 'export_annotated_features' => INPUT_PLUS() },
               'A->1' => { 'merge_annotated_features'  => INPUT_PLUS() }
            },
        },
        {   -logic_name  => 'export_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 10,
            -batch_size => 50,
            -parameters  => {
                cmd => 'export_annotated_features.pl --output_file #directory#/#file# --registry #reg_conf# --species #species# --min_id #min_id# --max_id #max_id#',
            },
        },
        {   -logic_name  => 'merge_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'cat #gff_files_from_batches# | xargs --max-args=50 cat >> #merged_gff#',
            },
            -flow_into   => {
               MAIN => 'gzip_annotated_features',
            },
        },
        {   -logic_name  => 'gzip_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #merged_gff#',
            },
            -flow_into   => {
               MAIN => 'mv_annotated_features_to_ftp',
            },
        },
        {   -logic_name  => 'mv_annotated_features_to_ftp',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #merged_gff#.gz #ftp_base_dir#/#species#',
            },
            -flow_into   => {
               MAIN => 'rm_annotated_features_temp_dir',
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
