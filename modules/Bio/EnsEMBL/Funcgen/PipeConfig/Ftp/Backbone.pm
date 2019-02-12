package Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::Backbone;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::Base';
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
               '1->A' => 'start_create_meta_files',
               'A->1' => 'backbone_fire_ftp_site_checker'
            },
        },
        {   -logic_name  => 'start_create_meta_files',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_ftp_site_checker',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_ftp_site_checker',
               'A->1' => 'backbone_ftp_pipeline_finished'
            },
        },
        {   -logic_name  => 'start_ftp_site_checker',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name => 'backbone_ftp_pipeline_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ]
}

1;
