package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportDataFilesAlignment;

use strict;
use warnings;
use Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportDataFilesBase;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportDataFilesBase';

sub pipeline_analyses {
    my $self = shift;

    my $dbfile_registry_path_parameter = $self->create_dbfile_registry_path_parameter;

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'fetch_parameters_from_db'
            },
        },
        {   -logic_name  => 'fetch_parameters_from_db',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::FetchParametersFromDb',
            -flow_into   => {
               MAIN => [
                  'export_data_files_bam',
                  'export_data_files_bigwig',
               ]
            },
        },
        {   -logic_name  => 'export_data_files_bam',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => '
                export_data_files.pl  \
                  --destination_root_path #ftp_base_dir#/#species#/Alignments  \
                  --file_type bam  \
                  --assembly #assembly# \
                ' . $dbfile_registry_path_parameter . ' \
                  --registry #reg_conf#  \
                  --die_if_source_files_missing 0 \
                  --data_freeze_date #data_freeze_date# \
                  --species #species#
                ',
            },
        },
        {   -logic_name  => 'export_data_files_bigwig',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => '
                export_data_files.pl  \
                  --destination_root_path #ftp_base_dir#/#species#/Alignments  \
                  --file_type bigwig  \
                  --assembly #assembly# \
                ' . $dbfile_registry_path_parameter . ' \
                  --registry #reg_conf#  \
                  --die_if_source_files_missing 0 \
                  --data_freeze_date #data_freeze_date# \
                  --species #species#
                ',
            },
        },
    ]
}

1;

