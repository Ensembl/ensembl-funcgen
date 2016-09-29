package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportQualityCheckValues;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Hash::Util qw( lock_hash );

sub pipeline_analyses {
    my $self = shift;

    my $data_freeze_date = $self->o('data_freeze_date');
    
    my $ftp_layout_configuration = {
      qc_chance_dir                       => '#ftp_base_dir#/#species#/QualityChecks',
      qc_mapped_reads_dir                 => '#ftp_base_dir#/#species#/QualityChecks',
      qc_phantom_peaks_dir                => '#ftp_base_dir#/#species#/QualityChecks',
      qc_proportion_of_reads_in_peaks_dir => '#ftp_base_dir#/#species#/QualityChecks',
      
      qc_chance_file_base_name                        => "#species#.#assembly#.chance.quality_check.${data_freeze_date}.json",
      qc_mapped_reads_base_name                       => "#species#.#assembly#.mapped_reads.quality_check.${data_freeze_date}.json",
      qc_phantom_peaks_file_base_name                 => "#species#.#assembly#.phantom_peaks.quality_check.${data_freeze_date}.json",
      qc_proportion_of_reads_in_peaks_file_base_name  => "#species#.#assembly#.proportion_of_reads_in_peaks.quality_check.${data_freeze_date}.json",
    };

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'job_factory_quality_checks'
            },
        },
        {   -logic_name  => 'job_factory_quality_checks',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryQualityChecks',
            -flow_into   => {
               MAIN => [
                'export_qc_chance',
                'export_qc_mapped_reads',
                'export_qc_phantom_peaks',
                'export_qc_proportion_of_reads_in_peaks'
              ]
            },
        },
        {   -logic_name  => 'export_qc_chance',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_qc_chance.pl --output_file ' 
                  . $ftp_layout_configuration->{qc_chance_dir} . '/' . $ftp_layout_configuration->{qc_chance_file_base_name} 
                  . ' --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_chance'
            },
        },
        {   -logic_name  => 'gzip_qc_chance',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip ' . $ftp_layout_configuration->{qc_chance_dir} . '/' . $ftp_layout_configuration->{qc_chance_file_base_name},
            },
        },
        {   -logic_name  => 'export_qc_mapped_reads',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_qc_mapped_reads.pl --output_file ' 
                  . $ftp_layout_configuration->{qc_mapped_reads_dir} . '/' . $ftp_layout_configuration->{qc_mapped_reads_base_name} 
                  . ' --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_mapped_reads'
            },
        },
        {   -logic_name  => 'gzip_qc_mapped_reads',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip ' . $ftp_layout_configuration->{qc_mapped_reads_dir} . '/' . $ftp_layout_configuration->{qc_mapped_reads_base_name},
            },
        },
        {   -logic_name  => 'export_qc_phantom_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_qc_phantom_peaks.pl --output_file ' 
                  . $ftp_layout_configuration->{qc_phantom_peaks_dir} . '/' . $ftp_layout_configuration->{qc_phantom_peaks_file_base_name} 
                  . ' --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_phantom_peaks'
            },
        },
        {   -logic_name  => 'gzip_qc_phantom_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip ' . $ftp_layout_configuration->{qc_phantom_peaks_dir} . '/' . $ftp_layout_configuration->{qc_phantom_peaks_file_base_name},
            },
        },
        {   -logic_name  => 'export_qc_proportion_of_reads_in_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_qc_proportion_of_reads_in_peaks.pl --output_file ' 
                  . $ftp_layout_configuration->{qc_proportion_of_reads_in_peaks_dir} . '/' . $ftp_layout_configuration->{qc_proportion_of_reads_in_peaks_file_base_name} 
                  . ' --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_proportion_of_reads_in_peaks'
            },
        },
        {   -logic_name  => 'gzip_qc_proportion_of_reads_in_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip ' . $ftp_layout_configuration->{qc_proportion_of_reads_in_peaks_dir} . '/' . $ftp_layout_configuration->{qc_proportion_of_reads_in_peaks_file_base_name},
            },
        },
    ]
}

1;
