package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportMotifFeatures;

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
               'MAIN' => 'export_motif_features'
            },
        },
        {   -logic_name  => 'export_motif_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_motif_features.pl --ftp_base_dir #ftp_base_dir#/#species# --registry #reg_conf# --species #species#',

            },
        },
    ]
}

1;
