package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Backbone;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'backbone_fire_exports',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_export',
               'A->1' => 'backbone_fire_createmetafiles'
            },
        },
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_createmetafiles',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'create_manifest_file',
               'A->1' => 'backbone_ftp_pipeline_finished'
            },
        },
        {   -logic_name  => 'create_manifest_file',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name => 'backbone_ftp_pipeline_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ]
}

1;
