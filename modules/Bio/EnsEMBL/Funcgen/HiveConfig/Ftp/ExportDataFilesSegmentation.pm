package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportDataFilesSegmentation;

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
                  'export_data_files_segmentation'
               ]
            },
        },
        {   -logic_name  => 'export_data_files_segmentation',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => '
                export_segmentation_files.pl            \
                  --registry #reg_conf#                 \
                  --species #species#                   \
                  --assembly #assembly#                 \
                  --data_freeze_date #data_freeze_date# \
                ' . $dbfile_registry_path_parameter . ' \
                  --destination_root_path #ftp_base_dir#/#species#/Segmentation  \
                ',
            },
        },
    ]
}

1;

