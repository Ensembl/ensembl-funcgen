package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportQualityCheckValues;

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
               'MAIN' => [
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
                cmd => 'export_qc_chance.pl --output_file #ftp_base_dir#/#species#/qc_chance.json --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_chance'
            },
        },
        {   -logic_name  => 'gzip_qc_chance',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #ftp_base_dir#/#species#/qc_chance.json',
            },
        },
        {   -logic_name  => 'export_qc_mapped_reads',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_qc_mapped_reads.pl --output_file #ftp_base_dir#/#species#/qc_mapped_reads.json --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_mapped_reads'
            },
        },
        {   -logic_name  => 'gzip_qc_mapped_reads',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #ftp_base_dir#/#species#/qc_mapped_reads.json',
            },
        },
        {   -logic_name  => 'export_qc_phantom_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_qc_phantom_peaks.pl --output_file #ftp_base_dir#/#species#/qc_phantom_peaks.json --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_phantom_peaks'
            },
        },
        {   -logic_name  => 'gzip_qc_phantom_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #ftp_base_dir#/#species#/qc_phantom_peaks.json',
            },
        },
        {   -logic_name  => 'export_qc_proportion_of_reads_in_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_qc_proportion_of_reads_in_peaks.pl --output_file #ftp_base_dir#/#species#/qc_proportion_of_reads_in_peaks.json --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_qc_proportion_of_reads_in_peaks'
            },
        },
        {   -logic_name  => 'gzip_qc_proportion_of_reads_in_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #ftp_base_dir#/#species#/qc_proportion_of_reads_in_peaks.json',
            },
        },
    ]
}

1;
