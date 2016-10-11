package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportMotifFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    my $data_freeze_date = $self->o('data_freeze_date');
    
    my $ftp_layout_configuration = {
      motif_features_dir             => '#ftp_base_dir#/#species#',
      motif_features_file_base_name  => "#species#.#assembly#.motiffeatures.${data_freeze_date}.gff",
    };

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'job_factory_motif_features'
            },
        },
        {   -logic_name  => 'job_factory_motif_features',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryQualityChecks',
            -flow_into   => {
               MAIN => 'export_motif_features'
            },
        },
        {   -logic_name  => 'export_motif_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_motif_features.pl --output_file '
                . $ftp_layout_configuration->{motif_features_dir} . '/' . $ftp_layout_configuration->{motif_features_file_base_name} 
                . ' --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_motif_features'
            },
        },
        {   -logic_name  => 'gzip_motif_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip ' . $ftp_layout_configuration->{motif_features_dir} . '/' . $ftp_layout_configuration->{motif_features_file_base_name} ,
            },
        },
    ]
}

1;
